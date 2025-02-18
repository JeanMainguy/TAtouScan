import pytest
from pathlib import Path
from tatouscan.annotation import find_ta_hits
from tatouscan.models import TaHit
from unittest import mock
from typing import Generator, List


class IterableMock(mock.Mock):
    def __init__(self, hits:List[mock.Mock]):
        super().__init__()
        self.hits = hits

    def __iter__(self):
        # Return an iterator that yields hits
        for hit in self.hits:
            yield hit

@pytest.fixture
def faa_file(tmp_path: Path) -> Path:
    """
    Fixture for a protein FASTA file.
    """
    faa_file: Path = tmp_path / "test.faa"
    faa_file.write_text(""">protein1
    MTEITAAMVKELRESTGAGMMDCKNALSETQHEFSQVLKAKLAEQAERYDDMAAAMKAVTEQGHELSNEER
    """)
    return faa_file


@pytest.fixture
def hmm_db(tmp_path: Path) -> Path:
    """
    Fixture for an HMM database file.
    """
    hmm_db: Path = tmp_path / "test.hmm"
    hmm_db.write_text("""HMMER3/f [3.1b2 | March 2015]
    NAME  TestHMM
    ACC   TEST00001
    DESC  Test HMM
    LENG  50
    ALPH  amino
    """)
    return hmm_db

@pytest.fixture
def mock_hits() -> IterableMock:

    mock_hit = mock.Mock()
    mock_hit.name = b"protein1"
    mock_hit.score = 100.0
    mock_hit.evalue = 0.001
    
    mock_query = mock.Mock()
    mock_query.name = b"TestHMM"
    
    mock_hits = IterableMock([mock_hit])
    mock_hits.query = mock_query
        
    return mock_hits


@pytest.fixture
def mock_hmmsearch(mock_hits:IterableMock) -> Generator[None, None, None]:
    """
    Fixture for mocking pyhmmer functionalities.
    """

    with mock.patch('pyhmmer.hmmsearch') as mock_hmmsearch:
        

        mock_hmmsearch.return_value = [mock_hits]

        yield mock_hmmsearch


def test_find_ta_hits_valid_input(faa_file: Path, hmm_db: Path, mock_hmmsearch:mock.Mock) -> None:
    """
    Test find_ta_hits with valid inputs and one hit, using a mock for hmmsearch.
    """

    
    result = find_ta_hits(faa_file, hmm_db)
    assert result is not None
    assert len(result) == 1
    assert "protein1" in result
    assert len(result["protein1"]) == 1
    assert isinstance(result["protein1"][0], TaHit)
    assert result["protein1"][0].protein_id == "protein1"
    assert result["protein1"][0].ta_name == "TestHMM"
    assert result["protein1"][0].score == 100.0
    assert result["protein1"][0].evalue == 0.001

def test_find_ta_hits_no_sequences(tmp_path:Path, hmm_db: Path
) -> None:
    """
    Test find_ta_hits when no sequences are in the FASTA file.
    """
        
    empty_faa_file: Path = tmp_path / "empty.faa"
    empty_faa_file.write_text("")

    with pytest.raises(ValueError):
        find_ta_hits(empty_faa_file, hmm_db)


def test_find_ta_hits_no_hits(
    faa_file: Path, hmm_db: Path, mock_hmmsearch: mock.Mock
) -> None:
    """
    Test find_ta_hits when no hits are found.
    """
    # Mock no hits
    mock_hmmsearch.return_value = []

    result = find_ta_hits(faa_file, hmm_db)
    assert result == {}
