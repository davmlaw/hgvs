"""Repeat (e.g. ``AGG[2]``) support: expanding parsed repeats into plain edits, and
detecting repeats in ref/alt sequences for pretty printing."""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

import hgvs
import hgvs.edit
import hgvs.location
import hgvs.posedit
import hgvs.sequencevariant
from hgvs.exceptions import HGVSInvalidVariantError, HGVSUnsupportedOperationError

if TYPE_CHECKING:
    from hgvs.pretty.models import VariantCoords


def expand_repeat(var, hdp):
    """Rewrite a repeat variant as the equivalent dup / ins / del / identity variant.

    HGVS repeat notation gives the total copy number of a repeat unit in the
    variant allele, e.g. ``NM_000059.4:c.938AGG[2]`` means "two copies of AGG
    starting at c.938". The reference is consulted to count the copies present
    at that position; the difference is then written as a plain edit placed at
    the 3' end of the repeat, following the HGVS 3' rule:

    * more copies than the reference: ``dup`` of the trailing copies when the
      extra copies fit inside the reference repeat, otherwise ``ins`` after it
    * fewer copies: ``del`` of the trailing copies
    * the same number: identity (``=``)

    Accepted forms: a start position with a unit (``c.938AGG[2]``), an interval
    with a unit (``c.6955_6993CAG[26]``), or an interval without a unit
    (``c.220_222[1]``), where the interval is one or more copies of the unit.

    Variants that are not repeats are returned unchanged. Uncertain copy
    numbers (``GCA[10_36]``, ``CGG[(90_?)]``) and intronic positions raise
    :class:`HGVSUnsupportedOperationError`. A unit that is not found in the
    reference at the stated position raises :class:`HGVSInvalidVariantError`.

    :param var: a :class:`hgvs.sequencevariant.SequenceVariant` of type c, g, m or n
    :param hdp: a data provider (for ``get_seq`` and, for c. variants, ``get_tx_identity_info``)
    :returns: a new :class:`hgvs.sequencevariant.SequenceVariant`
    """
    if var.posedit is None or not isinstance(var.posedit.edit, hgvs.edit.Repeat):
        return var

    edit = var.posedit.edit
    if edit.uncertain or edit.max is None or edit.min != edit.max:
        msg = f"Cannot expand repeat with uncertain copy number: {var}"
        raise HGVSUnsupportedOperationError(msg)
    if var.type not in "cgmn":
        msg = f"Repeat expansion is only supported for c, g, m and n variants: {var}"
        raise HGVSUnsupportedOperationError(msg)

    pos = var.posedit.pos
    if (
        pos is None
        or pos.uncertain
        or not _is_definite_position(pos.start)
        or not _is_definite_position(pos.end)
    ):
        msg = f"Cannot expand repeat with uncertain position: {var}"
        raise HGVSUnsupportedOperationError(msg)

    # Work in sequence (n. or g.) coordinates; convert c. positions via the CDS offset.
    mapper = None
    if var.type == "c":
        from hgvs.alignmentmapper import AlignmentMapper  # noqa: PLC0415 - import cycle

        mapper = AlignmentMapper(hdp, var.ac, var.ac, "transcript")
        pos = mapper.c_to_n(pos)
    if getattr(pos.start, "offset", 0) != 0 or getattr(pos.end, "offset", 0) != 0:
        msg = f"Repeat expansion of intronic variants is not supported: {var}"
        raise HGVSUnsupportedOperationError(msg)

    start_i = pos.start.base - 1  # 0-based, inclusive
    end_i = pos.end.base  # 0-based, exclusive

    unit = edit.ref.upper() if edit.ref else hdp.get_seq(var.ac, start_i, end_i).upper()
    unit_len = len(unit)
    span = end_i - start_i
    if span > 1 and span % unit_len != 0:
        msg = f"Repeat interval length ({span}) is not a multiple of the unit length ({unit_len}): {var}"
        raise HGVSInvalidVariantError(msg)

    ref_copies = _count_reference_copies(hdp, var.ac, start_i, unit)
    if ref_copies == 0:
        found = hdp.get_seq(var.ac, start_i, start_i + unit_len)
        msg = f"Repeat unit {unit} not found at position {pos.start} of {var.ac} (found {found}): {var}"
        raise HGVSInvalidVariantError(msg)
    if span > 1 and span // unit_len > ref_copies:
        msg = (
            f"Repeat interval spans {span // unit_len} copies of {unit} but the reference "
            f"has {ref_copies} at position {pos.start} of {var.ac}: {var}"
        )
        raise HGVSInvalidVariantError(msg)

    alt_copies = edit.min
    repeat_end = (
        pos.start.base + ref_copies * unit_len - 1
    )  # 1-based last base of the reference repeat
    if alt_copies == ref_copies:
        new_start, new_end = pos.start.base, repeat_end
        new_edit = hgvs.edit.NARefAlt(ref=unit * ref_copies, alt=unit * ref_copies)
    elif alt_copies > ref_copies:
        extra = alt_copies - ref_copies
        if extra <= ref_copies:
            new_start, new_end = repeat_end - extra * unit_len + 1, repeat_end
            new_edit = hgvs.edit.Dup(ref=unit * extra)
        else:
            new_start, new_end = repeat_end, repeat_end + 1
            new_edit = hgvs.edit.NARefAlt(ref=None, alt=unit * extra)
    else:
        removed = ref_copies - alt_copies
        new_start, new_end = repeat_end - removed * unit_len + 1, repeat_end
        new_edit = hgvs.edit.NARefAlt(ref=unit * removed, alt=None)

    start_pos = copy.deepcopy(pos.start)
    end_pos = copy.deepcopy(pos.end)
    start_pos.base = new_start
    end_pos.base = new_end
    new_pos = type(pos)(start=start_pos, end=end_pos, uncertain=False)
    if mapper is not None:
        new_pos = mapper.n_to_c(new_pos)

    return hgvs.sequencevariant.SequenceVariant(
        ac=var.ac,
        type=var.type,
        posedit=hgvs.posedit.PosEdit(pos=new_pos, edit=new_edit),
        gene=var.gene,
    )


def _is_definite_position(p):
    return (
        p is not None
        and not isinstance(p, hgvs.location.Interval)
        and getattr(p, "base", None) is not None
        and not p.uncertain
    )


def _count_reference_copies(hdp, ac, start_i, unit, chunk_units=64):
    """Count consecutive copies of `unit` in sequence `ac` beginning at 0-based `start_i`."""
    unit_len = len(unit)
    copies = 0
    pos = start_i
    while True:
        chunk = hdp.get_seq(ac, pos, pos + unit_len * chunk_units).upper()
        n = 0
        while (n + 1) * unit_len <= len(chunk) and chunk.startswith(unit, n * unit_len):
            n += 1
        copies += n
        if n < chunk_units:
            return copies
        pos += unit_len * chunk_units


@dataclass(eq=True, repr=True, frozen=True, order=True)
class RepeatUnit:
    repeat_count: int
    repeat_unit: str
    block_size: int
    block: str


class RepeatAnalyser:
    def __init__(self, fs: VariantCoords, reverse: bool = False) -> None:
        self.is_repeat = False
        self.repeat_units_ref = None
        self.repeat_units_alt = None
        self.ref_str = fs.ref
        self.alt_str = fs.alt
        self.reverse = reverse

        # protect from too large sequences getting analyzed.
        if len(self.ref_str) > hgvs.global_config.repeats.max_repeat_length:
            self.repeat_units_ref = []
            self.repeat_units_alt = []
            return

        self.repeat_units_ref = detect_repetitive_block_lengths(self.ref_str, reverse=self.reverse)
        self.repeat_units_alt = detect_repetitive_block_lengths(self.alt_str, reverse=self.reverse)

        if len(self.repeat_units_ref) == 0 and len(self.repeat_units_alt) == 0:
            return

        # check longest repeat blocks:
        # we only look at ref to determine if there are repeats
        # If ref has no repeat, we don't call this a repeat variant, even if alt would have a repetitive unit
        longest_r_unit = self._get_longest_repeat_unit(self.repeat_units_ref)
        if longest_r_unit is None:
            return

        # filter our too fragmented results
        expected_size = len(self.ref_str) / 3
        if longest_r_unit.block_size < expected_size:
            return

        if longest_r_unit.repeat_unit not in self.alt_str:
            return

        self.repeat_units_alt = detect_repetitive_block_lengths(
            self.alt_str, longest_ref_unit=longest_r_unit, reverse=self.reverse
        )

        self.is_repeat = True

        ref_repeat = get_repeat_str(
            self.ref_str, self.alt_str, self.repeat_units_ref, self.repeat_units_alt, self.reverse
        )
        alt_repeat = get_repeat_str(
            self.alt_str, self.ref_str, self.repeat_units_alt, self.repeat_units_ref, self.reverse
        )
        self.ref_str = ref_repeat
        self.alt_str = alt_repeat

    def __repr__(self):
        return f"{self.ref_str}>{self.alt_str}"

    def _get_longest_repeat_unit(self, repeat_units: list[RepeatUnit]) -> RepeatUnit:
        lru = None
        for ru in repeat_units:
            if not lru:
                lru = ru
                continue

            if ru.block_size > len(lru.block):
                lru = ru

        return lru


def detect_repetitive_block_lengths(
    sequence: str, longest_ref_unit: RepeatUnit | None = None, reverse: bool = False
) -> list[RepeatUnit]:
    """Detects the length of repetitive blocks in a string, with an option to search from left to right or reverse.

    In reverse mode, it creates the largest possible blocks of the smallest possible units.
    """
    result: list[RepeatUnit] = []
    seq_len = len(sequence)

    if longest_ref_unit is not None:
        # look for full containment of the longest ref repeat
        # this is so we can detect [2]>[1] (or really any repeat length variation to just 1)
        pattern = f"({re.escape(longest_ref_unit.repeat_unit)})+"
        match = re.search(pattern, sequence)

        if match:
            repeat_count = len(sequence) // len(longest_ref_unit.repeat_unit)
            block = repeat_count * longest_ref_unit.repeat_unit
            ru = RepeatUnit(repeat_count, longest_ref_unit.repeat_unit, len(block), block)
            result.append(ru)
            shuffleable_bases = sequence[repeat_count * len(ru.repeat_unit) :]
            rus = detect_repetitive_block_lengths(shuffleable_bases, reverse=reverse)
            result.extend(rus)
            return result

    if reverse:
        i = seq_len  # Start from the end of the sequence
        while i > 0:
            matched = False
            # Iterate over block sizes from smallest to largest
            for block_size in range(1, seq_len // 2 + 1):
                if i - block_size < 0:
                    continue  # Not enough characters to form a repeat unit

                # Extract the potential repeat unit ending at position i
                repeat_unit = sequence[i - block_size : i]

                # Calculate the maximum possible number of repeats for this block size
                max_possible_repeats = i // block_size

                # Initialize repeat_count to the maximum possible and decrease
                for repeat_count in range(max_possible_repeats, 1, -1):
                    start_index = i - repeat_count * block_size
                    if start_index < 0:
                        continue  # Not enough characters to form the repetitive block

                    # Extract the substring that could be the repetitive block
                    substr = sequence[start_index:i]

                    # Build the regex pattern for the current repeat unit and count
                    pattern = rf"({re.escape(repeat_unit)})" + r"{" + f"{repeat_count}" + r"}"

                    # Check if the substring matches the pattern
                    if re.fullmatch(pattern, substr):
                        repetitive_block = substr
                        ru = RepeatUnit(
                            block=repetitive_block,
                            block_size=len(repetitive_block),
                            repeat_unit=repeat_unit,
                            repeat_count=repeat_count,
                        )
                        result.append(ru)
                        # Move the index `i` backward by the length of the repetitive block
                        i -= len(repetitive_block)
                        matched = True
                        break  # Found the largest block for this unit size

                if matched:
                    break  # Proceed to the next position `i`

            if not matched:
                # No repeat found, remove one character
                ru = RepeatUnit(
                    block=sequence[i - 1], block_size=1, repeat_unit=sequence[i - 1], repeat_count=1
                )
                result.append(ru)
                i -= 1  # Move back by one character if no match is found

    else:
        i = 0  # Start from the beginning of the sequence
        while i < seq_len:
            matched = False
            for block_size in range(1, seq_len // 2 + 1):
                # Build the regex pattern for the current block size
                pattern = rf"(.{{{block_size}}})\1+"
                match = re.match(pattern, sequence[i:])

                if match:
                    repetitive_block = match.group()  # The full repeating pattern
                    repeated_unit = match.group(1)  # The repeating unit
                    repetition_count = len(repetitive_block) // len(repeated_unit)

                    # Add the repetitive block and its details to the result
                    ru = RepeatUnit(
                        block=repetitive_block,
                        block_size=len(repetitive_block),
                        repeat_unit=repeated_unit,
                        repeat_count=repetition_count,
                    )
                    result.append(ru)

                    # Move the index `i` forward by the length of the repetitive block
                    i += len(repetitive_block)
                    matched = True
                    break  # Found a match, break the loop for current block size

            if not matched:
                i += 1  # Move forward by one character if no match is found

    return result


def get_repeat_str(
    sequence: str,
    other_seq: str,
    primary_repeat_unit: list[RepeatUnit],
    other_repeat_unit: list[RepeatUnit],
    reverse: bool = False,
) -> str:
    if len(primary_repeat_unit) == 0 and len(other_repeat_unit) == 0:
        return None

    if (
        len(primary_repeat_unit) == 0
        and len(other_repeat_unit) == 1
        and sequence == other_repeat_unit[0].repeat_unit
    ):
        return f"{sequence}[1]"
    if (
        len(primary_repeat_unit) == 0
        and len(other_repeat_unit) == 1
        and sequence != other_repeat_unit[0].repeat_unit
    ):
        return None

    if len(primary_repeat_unit) > 0 and len(other_repeat_unit) > 0:
        return_str = assemble_repeat_string(sequence, primary_repeat_unit, reverse=reverse)
        return return_str

    if len(other_repeat_unit) == 0 and len(other_seq) > 0:
        return_str = assemble_repeat_string(sequence, primary_repeat_unit, reverse=reverse)
        if len(return_str) > 0:
            return return_str

    return None


def assemble_repeat_string(
    sequence: str, repeat_units: list[RepeatUnit], reverse: bool = False
) -> str:
    return_str = ""
    primary_repeat_unit = repeat_units.copy()
    seq = sequence
    if reverse:
        while len(seq) > 0:
            found_unit = None
            for ru in primary_repeat_unit:
                if seq.endswith(ru.block):
                    return_str = f"{ru.repeat_unit}[{ru.repeat_count}]" + return_str
                    seq = seq[: -len(ru.block)]
                    found_unit = ru
                    break
            if not found_unit:
                # remove one character if no matching repeat unit is found
                count = 1
                seq_char = seq[0]
                seq = seq[:-1]

                # count consecutive repeating chars  from the end
                while seq and seq[-1] == seq_char:
                    count += 1
                    seq = seq[:-1]

                return_str = f"{seq_char}[{count}]" + return_str

    else:  # forward direction
        while len(seq) > 0:
            found_unit = None
            for ru in primary_repeat_unit:
                if seq.startswith(ru.block):
                    return_str += f"{ru.repeat_unit}[{ru.repeat_count}]"
                    seq = seq[len(ru.block) :]
                    found_unit = ru
                    break
            if not found_unit:
                # remove one character if no matching repeat unit is found
                count = 1
                seq_char = seq[0]
                seq = seq[1:]

                # count consecutive repeating chars
                while seq and seq[0] == seq_char:
                    count += 1
                    seq = seq[1:]

                return_str += f"{seq_char}[{count}]"

    return return_str
