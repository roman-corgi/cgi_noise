import os
import csv
import pandas as pd
import chardet
from pathlib import Path

def open_folder(*folders):
    """Opens a directory and returns a dictionary of file paths keyed by filenames."""
    filenamedir = Path(os.getcwd()).parent
    folder = Path(filenamedir, *folders)
    return {file.name: file for file in folder.iterdir() if file.is_file()}

class loadCSVrow:
    """
    Class to load CSV data into a pandas DataFrame, stopping at rows containing any cell that starts with a '#'.

    Parameters
    ----------
    filename : str
        Path to the CSV file.
    
    Attributes
    ----------
    df : pd.DataFrame
        DataFrame containing CSV data, excluding comment rows and rows with in-line comments.
    comments : list
        List of comments extracted from the file.
    """
    
    def __init__(self, filename):
        self.fullfile = filename
        self.csvname = os.path.basename(filename)
        self.prefix = os.path.basename(filename).split("_")[0]
        self._encoding = self.detect_encoding()
        self.comments = []
        self.df = self.load_data()

    def detect_encoding(self):
        """Detect file encoding."""
        with open(self.fullfile, 'rb') as f:
            result = chardet.detect(f.read())
            return result['encoding']

    def load_data(self):
        """Loads data into a DataFrame, converting each cell to numeric if possible, while preserving text data."""
        rows = []
        is_comment_section = False

        with open(self.fullfile, 'r', encoding=self._encoding) as f:
            reader = csv.reader(f)
            for row in reader:
                # If a comment line is encountered, switch to comment mode
                if row and row[0].startswith("#"):
                    is_comment_section = True
                
                # Append all rows to comments if in comment section
                if is_comment_section:
                    self.comments.append(row)
                else:
                    # Check if any cell in the row is a comment (starts with '#')
                    if any(cell.startswith("#") for cell in row):
                        break  # Stop reading further if an inline comment is detected
                    rows.append(row)

        # Convert collected rows to a DataFrame (excluding header rows if they are comments)
        if rows:
            header = rows[0]  # Assuming first row is the header
            data = rows[1:]   # Data starts after header
            df = pd.DataFrame(data, columns=header, dtype=str)  # Load everything as strings initially

            # Convert each cell to numeric if possible
            for col in df.columns:
                df[col] = df[col].apply(self._convert_cell_to_numeric)

        else:
            df = pd.DataFrame()  # Return an empty DataFrame if no data rows are collected

        return df

    def _convert_cell_to_numeric(self, cell):
        """Helper function to convert a cell to a numeric type if possible."""
        try:
            # Try converting to a float first (covers both integers and decimals)
            return float(cell)
        except ValueError:
            # If conversion fails, return the cell as-is (keeps it as text)
            return cell
