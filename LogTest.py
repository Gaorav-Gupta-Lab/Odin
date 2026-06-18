# myapp.py
import logging

def main():
    log_file_name = "/mnt/hgfs/D_Drive/Testing/Test.log"
    logging.basicConfig(filename=log_file_name, level=logging.INFO)
    logging.info('Started')

    logging.info('Finished')

if __name__ == '__main__':
    main()