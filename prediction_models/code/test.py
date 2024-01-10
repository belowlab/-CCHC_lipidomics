import logging
# f string is not fully compatible with logging, so use %s for string formmatting
logging.root.handlers = []
logging.basicConfig(level=logging.DEBUG, handlers=[logging.FileHandler(filename='app.log', mode='w'), logging.StreamHandler()], format='%(name)s - %(levelname)s - %(message)s')
# logging.basicConfig(level=logging.DEBUG, filename='app.log', filemode='a', format='%(name)s - %(levelname)s - %(message)s')
logging.debug('This is a debug message')
logging.info('This is an info message')
logging.warning('This is a warning message')
logging.error('This is an error message')
logging.critical('This is a critical message')


