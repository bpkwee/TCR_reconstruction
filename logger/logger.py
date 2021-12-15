import logging
import logging.handlers
import os


def init_logger(log_path_and_file='log_file_UNKNOWN_script.log', level_msg='DEBUG'):
    should_roll_over = os.path.isfile(log_path_and_file)

    if level_msg == 'DEBUG':
        level_msg = logging.DEBUG
    elif level_msg == 'INFO':
        level_msg = logging.INFO
    elif level_msg == 'WARNING':
        level_msg = logging.WARNING

    logging.basicConfig(filename=log_path_and_file,
                        level=level_msg,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S',
                        filemode='w')

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', '%m/%d/%Y %I:%M:%S')
    handler = logging.handlers.RotatingFileHandler(log_path_and_file,
                                                   mode='w',
                                                   backupCount=0)  # ensure older backups get deleted

    handler.setFormatter(formatter)

    log = logging.getLogger(log_path_and_file)
    log.addHandler(handler)

    if should_roll_over:  # log already exists, roll over!
        handler.doRollover()

    return log
