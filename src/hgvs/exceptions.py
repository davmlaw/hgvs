"""defined exceptions used by the hgvs package"""


class HGVSError(Exception):
    pass


class HGVSDataNotAvailableError(HGVSError):
    pass


class HGVSInternalError(HGVSError):
    pass


class HGVSInvalidIntervalError(HGVSError):
    pass


class HGVSInvalidVariantError(HGVSError):
    """Exception raised when variant is inconsistent or invalid"""


class HGVSExonJunctionSpanError(HGVSInvalidVariantError):
    """Exception raised when a transcript (c./n.) interval spans an exon junction: its ends
    lie in different exons, so the bases it covers are adjacent on the transcript but
    separated by an intron on the genome, and the variant has no single genomic equivalent"""


class HGVSNormalizationError(HGVSError):
    pass


class HGVSParseError(HGVSError):
    pass


class HGVSUnsupportedOperationError(HGVSError):
    pass


class HGVSUsageError(HGVSError):
    """Exception raised when client/caller has made an invalid request"""


class HGVSVerifyFailedError(HGVSError):
    """Exception raised when the cached hdp result is not consistent with latest remote result"""


# <LICENSE>
# Copyright 2018 HGVS Contributors (https://github.com/biocommons/hgvs)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
