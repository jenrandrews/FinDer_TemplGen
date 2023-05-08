DEFLOG = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'fileHandler': {
            'level': 'INFO',
            'formatter': 'fileFormatter',
            'class': 'logging.FileHandler',
            'mode': 'w',
        },
    },
    'formatters': {
        'fileFormatter': {
            'format': '%(asctime)s - %(levelname)s -  %(name)s - %(module)s - %(funcName)s - %(message)s',
            'datefmt': '%Y-%m-%d %H:%M:%S',
            'class': 'logging.Formatter'
        },
    },
    'loggers': {
        '': {
            'handlers': ['fileHandler'],
            'level': 'INFO',
        },
    }
}
