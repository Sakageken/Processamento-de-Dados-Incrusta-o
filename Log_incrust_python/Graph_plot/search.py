import numpy as np
from pathlib import Path
from scipy.fft import fft # Using scipy's fft

def check_readen_files(txt_filepath, filename_to_check):
    """
    This function serves to run through a txt file and search if 
    there is a specified text inside it, returns True if the file is there
    and false if its not
    """
    #Checks if a filename exists as a distinct line in the log text file.
    if Path(txt_filepath).is_file():
        pass
    else:
        #If the file does not exists, create all the folder necessary for its existence
        target_directory = Path.home()/"Downloads"
        target_directory.mkdir(parents=True, exist_ok=True)
        with open(txt_filepath, 'w') as f:
            f.write("")
        print("File .txt generated")

    try:
        #Open the txt file and check if the name provided are inside the txt file
        #line by line
        with open(txt_filepath, 'r') as f:
            for line in f:
                # .strip() removes leading/trailing whitespace, including the newline character
                if line.strip() == filename_to_check:
                    #If the name is there return True
                    return True
        # Return False if the loop finishes without finding a match
        return False 
    #If no file txt is found returns an error
    except FileNotFoundError:
        print(f"Error: The file '{txt_filepath}' was not found.")
        return False

def search_data(path_files: str):
    """
        This function searchs in the download home path for the Log_files.xlsx
        that exists and checks if they have already been read or not, 
        it returns a list with all the logs not read
    """
    try:
        #Searching for logs in the path_files folder
        search_folder = Path(Path.home() / path_files)
        #Path to the documents home path
        path_analised_docs = Path(Path.home()/"Documents")
        #Path to where all the already read logs are
        destination_logs = path_analised_docs / "Logs_analysed" / "Excel_files"
        destination = path_analised_docs / "Logs_analysed"
        #Path to the txt file where all tha names of the already read logs are
        path_txt_file = destination / "log_readen.txt"
        #Array with all the pathname of the logs in the download folder
        log_file = list(search_folder.glob("Data Log*.xlsx"))
        #Array to insert the files that will be read
        list_files_to_read = []
        #Run through the list of the files found and check if they 
        # have already been read
        for qtd in range(0, len(log_file)):
            #Check if destination exists
            if destination.is_dir():
                pass
            else:
                #If it did not then create all the folder for its existence
                destination.mkdir(parents=True, exist_ok=True)
            #Check if destination_logs exists
            if destination_logs.is_dir():
                pass
            else:
                #If it did not then create all the folder for its existence
                destination_logs.mkdir(parents=True, exist_ok=True)
            #Check for the name inside the txt file
            search = check_readen_files(f"{path_txt_file}", log_file[qtd].name)
            # If file already readen then move it to another folder
            if search:
                #Assign the old path and name to a variable
                file_to_be_moved = Path(log_file[qtd])
                file_name = log_file[qtd].name
                try:
                    #Tries to rename the file pathname so it can move it
                    file_to_be_moved.rename(destination_logs/file_name)
                except FileNotFoundError:
                    #If the file did not exist then print an error message
                    print("Error: The source file or destination folder does not exist.")
            else:
                #If the file was not in the txt file then add to a list 
                # so it can be read
                print("File not readen, should be readen and then moved to a new folder")
                list_files_to_read.append(log_file[qtd])

        #Return the list of files to be read, txt filepath, folder to store the read files
        return np.array([list_files_to_read]), path_txt_file, destination
    
    except FileNotFoundError:
        print("No file with the prefix 'Data Log' found in the folder /Users/lucasamorim/Downloads")
        return 

def write_on_text(txt_path: str, frase: str):
    """
    This function writes a sentece inside a txt file and break 
    the line by the end of the run, does not return anything
    """
    try:
        txt_path = Path(txt_path)
        with open(txt_path, 'a') as f:
            f.write(f"{frase}\n")
            print(f"Log: {frase} written on {txt_path} file")
    except FileNotFoundError:
        print(f"No txt file found in {txt_path}")
        

if __name__ == "__main__":
    print("File should not be run as main")