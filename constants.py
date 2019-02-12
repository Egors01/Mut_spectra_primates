class MutConstants():
    def __init__(self):
        self.chr_list = ['chr1', 'chr2', 'chr3', 'chr4',
             'chr5', 'chr6', 'chr7', 'chr8',
             'chr9', 'chr10', 'chr11', 'chr12',
             'chr13', 'chr14', 'chr15', 'chr16',
             'chr17', 'chr18', 'chr19', 'chr20',
             'chr21', 'chr22']
        self.name_to_colname_example_dict = {
        'Gorilla': 'gorilla',
        'Orang': 'orang',
        'Chimp':'chimp',
        'Bonobo':'bonobo',
        'Human':'human',
        'Gorilla_P': 'gorilla_P',
        'Orang_P': 'orang_P',
        'Chimp_P':'chimp_P',
        'Bonobo_P':'bonobo_P',
        'Human_P':'human_P'}

        self.pictures_pol_example = [
            [['Human_P','Chimp_P'],
            ['Human_P','Orang_P'],
            ['Chimp_P','Orang_P'],
            ['Bonobo_P','Orang_P'],
            ['Gorilla_P','Orang_P']],

            [['Human_P','Chimp_P'],
            ['Bonobo_P','Chimp_P'],
            ['Gorilla_P','Chimp_P'],
            ['Orang_P','Chimp_P']],

            [['Human_P','Gorilla_P'],
            ['Chimp_P','Gorilla_P'],
            ['Bonobo_P','Gorilla_P'],
            ['Orang_P','Gorilla_P']],

            [['Human_P','Bonobo_P'],
            ['Chimp_P','Bonobo_P'],
            ['Gorilla_P','Bonobo_P'],
            ['Orang_P','Bonobo_P']]
                                  ]
        self.pictures_div_example = [
            [['Human','Chimp'],
            ['Human','Orang'],
            ['Chimp','Orang'],
            ['Bonobo','Orang'],
            ['Gorilla','Orang']],

            [['Human','Chimp'],
            ['Bonobo','Chimp'],
            ['Gorilla','Chimp'],
            ['Orang','Chimp']],

            [['Human','Gorilla'],
            ['Chimp','Gorilla'],
            ['Bonobo','Gorilla'],
            ['Orang','Gorilla']],

            [['Human','Bonobo'],
            ['Chimp','Bonobo'],
            ['Gorilla','Bonobo'],
            ['Orang','Bonobo']],

             [['Human','Human_P'],
            ['Chimp','Chimp_P'],
            ['Bonobo','Bonobo_P'],
            ['Gorilla','Gorilla_P'],
            ['Orang','Orang_P']]
               ]

