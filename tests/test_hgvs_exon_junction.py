"""Transcript intervals that span an exon junction have no single genomic equivalent.

NM_000059.4:c.8331_8332del is the motivating case: c.8331 is the last base of one
BRCA2 exon and c.8332 the first base of the next, so mapping each end separately
yields a 6,870 bp genomic deletion that includes the intron. These tests use
transcripts in the offline test cache; explicit reference sequences avoid sequence
lookups that the cache does not hold.
"""

import pytest

import hgvs.assemblymapper
import hgvs.parsers
import hgvs.variantmapper
from hgvs.exceptions import HGVSExonJunctionSpanError, HGVSInvalidVariantError

# NM_000540.2 (RYR1, + strand on NC_000019.9): exons 10 and 11 (ord 9 and 10) meet at
# tx index 1087, i.e. n.1087 is the last base of one exon and n.1088 the first of the
# next. cds_start_i is 130, so n.1087 == c.957. Genomic: exon ends at g.38939151, the
# next starts at g.38939289.
PLUS_TX, PLUS_ALT = "NM_000540.2", "NC_000019.9"
# NM_000138.4 (FBN1, - strand on NC_000015.9): junction at tx index 1542 (n.1542 |
# n.1543; cds_start_i 395 so c.1147 | c.1148). Genomic: g.48812856 and g.48808559.
MINUS_TX, MINUS_ALT = "NM_000138.4", "NC_000015.9"


@pytest.fixture(scope="module")
def hp():
    return hgvs.parsers.Parser()


@pytest.fixture(scope="module")
def strict_vm(hdp):
    return hgvs.variantmapper.VariantMapper(hdp, prevalidation_level=None, replace_reference=False)


@pytest.fixture(scope="module")
def lenient_vm(hdp):
    return hgvs.variantmapper.VariantMapper(
        hdp, prevalidation_level=None, replace_reference=False, strict_exon_junction=False
    )


class TestAlignmentMapperJunctions:
    def test_junctions_recorded(self, strict_vm):
        mapper = strict_vm._fetch_AlignmentMapper(PLUS_TX, PLUS_ALT, "splign")
        assert len(mapper.tx_exon_junctions) == 105  # 106 exons
        assert 1087 in mapper.tx_exon_junctions
        assert mapper.tx_exon_junctions == sorted(mapper.tx_exon_junctions)

    def test_transcript_mapper_has_no_junctions(self, strict_vm):
        mapper = strict_vm._fetch_AlignmentMapper(PLUS_TX, PLUS_TX, "transcript")
        assert mapper.tx_exon_junctions == []

    @pytest.mark.parametrize(
        ("interval", "expected"),
        [
            ("1087_1088", [1087]),  # last base of exon | first base of the next
            ("1080_1095", [1087]),
            ("930_1088", [930, 1087]),  # two junctions (exon ord 9 is tx 930..1087)
            ("1086_1087", []),  # ends inside one exon
            ("1088_1089", []),
            ("1087", []),  # single base
            ("1087_1087+1", []),  # intronic offset: the end is in the intron, not the next exon
            ("1087+1_1088", []),
            ("(1080_1087)_(1088_1095)", []),  # uncertain ends are not judged
        ],
    )
    def test_n_exon_junctions_spanned(self, hp, strict_vm, interval, expected):
        mapper = strict_vm._fetch_AlignmentMapper(PLUS_TX, PLUS_ALT, "splign")
        assert mapper.n_exon_junctions_spanned(hp.parse_n_interval(interval)) == expected

    def test_minus_strand_uses_transcript_coordinates(self, hp, strict_vm):
        mapper = strict_vm._fetch_AlignmentMapper(MINUS_TX, MINUS_ALT, "splign")
        assert mapper.n_exon_junctions_spanned(hp.parse_n_interval("1542_1543")) == [1542]
        assert mapper.n_exon_junctions_spanned(hp.parse_n_interval("1541_1542")) == []


class TestVariantMapperStrict:
    @pytest.mark.parametrize(
        "hgvs_str",
        [
            "NM_000540.2:n.1087_1088delNN",
            "NM_000540.2:n.1087_1088delNNinsA",
            "NM_000540.2:n.1087_1088dupNN",
            "NM_000540.2:n.1087_1088invNN",
            "NM_000540.2:n.1087_1088insATGCGG",
            "NM_000540.2:n.1080_1095delNNNNNNNNNNNNNNNN",
        ],
    )
    def test_n_to_g_raises(self, hp, strict_vm, hgvs_str):
        with pytest.raises(HGVSExonJunctionSpanError, match=r"between n\.1087 and n\.1088"):
            strict_vm.n_to_g(hp.parse(hgvs_str), PLUS_ALT)

    def test_c_to_g_raises(self, hp, strict_vm):
        with pytest.raises(HGVSExonJunctionSpanError, match=r"NM_000540\.2:c\.957_958del"):
            strict_vm.c_to_g(hp.parse("NM_000540.2:c.957_958delNN"), PLUS_ALT)

    def test_c_to_g_raises_minus_strand(self, hp, strict_vm):
        with pytest.raises(HGVSExonJunctionSpanError, match=r"between n\.1542 and n\.1543"):
            strict_vm.c_to_g(hp.parse("NM_000138.4:c.1147_1148delNN"), MINUS_ALT)

    def test_error_is_an_invalid_variant_error(self, hp, strict_vm):
        with pytest.raises(HGVSInvalidVariantError):
            strict_vm.n_to_g(hp.parse("NM_000540.2:n.1087_1088delNN"), PLUS_ALT)

    @pytest.mark.parametrize(
        ("hgvs_str", "expected"),
        [
            # within one exon
            ("NM_000540.2:n.1086_1087delNN", "NC_000019.9:g.38939150_38939151del"),
            ("NM_000540.2:n.1088_1089delNN", "NC_000019.9:g.38939289_38939290del"),
            # intronic offsets describe genomic positions directly
            ("NM_000540.2:n.1087_1087+1delNN", "NC_000019.9:g.38939151_38939152del"),
            # substitutions and identity are not affected
            ("NM_000540.2:n.1087N>A", "NC_000019.9:g.38939151N>A"),
            ("NM_000540.2:n.1087_1088NN=", "NC_000019.9:g.38939151_38939289="),
        ],
    )
    def test_n_to_g_still_maps(self, hp, strict_vm, hgvs_str, expected):
        assert str(strict_vm.n_to_g(hp.parse(hgvs_str), PLUS_ALT)) == expected


class TestVariantMapperLenient:
    """strict_exon_junction=False keeps the historical behaviour: the genomic span
    between the two mapped ends, intron included (as VariantValidator also returns)."""

    def test_n_to_g_spans_intron(self, hp, lenient_vm):
        var_g = lenient_vm.n_to_g(hp.parse("NM_000540.2:n.1087_1088delNN"), PLUS_ALT)
        assert str(var_g) == "NC_000019.9:g.38939151_38939289del"

    def test_c_to_g_spans_intron_minus_strand(self, hp, lenient_vm):
        var_g = lenient_vm.c_to_g(hp.parse("NM_000138.4:c.1147_1148delNN"), MINUS_ALT)
        assert str(var_g) == "NC_000015.9:g.48808559_48812856del"

    def test_ins_spans_intron(self, hp, lenient_vm):
        var_g = lenient_vm.n_to_g(hp.parse("NM_000540.2:n.1087_1088insATGCGG"), PLUS_ALT)
        # historical output: the intron between the two exonic ends is replaced
        assert str(var_g) == "NC_000019.9:g.38939152_38939288delinsATGCGG"


def test_assemblymapper_passes_option_through(hdp):
    am = hgvs.assemblymapper.AssemblyMapper(
        hdp,
        assembly_name="GRCh37",
        normalize=False,
        prevalidation_level=None,
        replace_reference=False,
        strict_exon_junction=False,
    )
    assert am.strict_exon_junction is False
    am = hgvs.assemblymapper.AssemblyMapper(
        hdp, assembly_name="GRCh37", normalize=False, prevalidation_level=None
    )
    assert am.strict_exon_junction is True
