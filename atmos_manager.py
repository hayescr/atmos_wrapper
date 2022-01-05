import os
from abund_utils import atomic_sym_to_num


class AtmosManager:
    '''
    This class is provides a python interface for a\MARCS model atmosphere grid
        interpolator wrapper. It can also add MOOG formated chemical abundances
        to interpolated MARCS atmosphere files.
    '''

    wrapper_path = os.path.dirname(os.path.realpath(__file__))

    def __init__(self, teff, logg, feh, vmicro=None, cfe=0., alphafe=0.):
        '''
        On initialization the AtmosManager will set the input stellar
            parameters.

        teff : float, effective temperature
        logg : float, surface gravity
        feh : float, metallicity
        vmicro : float, microturbulent velocity
        cfe : float, carbon to iron ratio
        alphafe : float, alpha to iron ratio
        '''
        self.teff = teff
        self.logg = logg
        self.feh = feh
        if vmicro is None:
            self.vmicro = self.vmicro_calc()
        else:
            self.vmicro = vmicro

        self.cfe = cfe
        self.alphafe = alphafe

    def vmicro_calc(self):
        '''
        Custom vmicro calculation from BACCHUS.  Used if no vmicro is provided.

        returns

        vmicro : float, microturbulent velocity
        '''
        if self.logg < 3.5:
            if self.teff < 5250.:
                vmicro = 1.15 + 2e-4 * (self.teff - 5500.) + 3.95e-7 * (
                    self.teff - 5500.)**2. - 0.13 * (
                    self.logg - 4.0) + 0.13 * (self.logg - 4.0) ** 2.
            else:
                vmicro = 1.15 + 2e-4 * (5250. - 5500.) + 3.95e-7 * (
                    5250. - 5500.)**2. - 0.13 * (
                    self.logg - 4.0) + 0.13 * (self.logg - 4.0) ** 2.
        else:
            vmicro = 0.94 + 2.2E-5 * (self.teff - 5500) - 0.5E-7 * (
                self.teff - 5500)**2 - 0.1 * (
                self.logg - 4.0) + 0.04 * (
                self.logg - 4.0)**2 - 0.37 * self.feh - 0.07 * self.feh**2

        return vmicro

    def interp_atmos(self, star=None, file_format='ts', path=None):
        '''
        Formats the input for the wrapper of the MARCS model atmosphere
            interpolater.

        star : string, input star name (optional)
        file_format : string, indicates whether the interoplated atmospheres
            are to be in MOOG or turbospectrum formats.
        path : string, path to the output interpolated atmosphere

        returns

        filename : string, the output interpolated atmosphere filename
        '''
        if file_format.lower() in ['moog', 'm']:
            format_code = 'MOOG'
        elif file_format.lower() in ['turbospectrum',
                                     'turbo', 'turbospec', 'ts']:
            format_code = 'TurboSpectrum'

        if star is None:
            self.filename = (f'{self.teff:.0f}g{self.logg:.2f}m1.0z'
                             f'{self.feh:.2f}.int')
        else:
            self.filename = (f'{self.teff:.0f}g{self.logg:.2f}m1.0z'
                             f'{self.feh:.2f}_{star}.int')

        atmos_param_arg = (f'{self.teff:.0f} {self.logg:.2f} {self.feh:.2f} '
                           f'{self.vmicro:.2f} {self.cfe:.2f} '
                           f'{self.alphafe:.2f} {format_code} {self.filename}')

        cwd = os.getcwd()
        os.chdir(AtmosManager.wrapper_path)
        os.system(f'./load_atmos_param.com {atmos_param_arg}')

        if path is None:
            path = cwd

        if os.path.isfile(f'./{self.filename}'):
            os.rename(f'./{self.filename}', f'{path}/{self.filename}')
            self.file_loc = f'{path}/{self.filename}'
        else:
            print('Interpolation failed')

        os.chdir(cwd)

        return filename

    def moog_abund(self, abund_dict):
        '''
        Takes input dictionary of log epsilon abundances (e.g., {'C' : 8.66}),
            and appends the necessary chemical abundance information in MOOG
            format to the interpolated model atmosphere.

        abund_dict : dictionary, contains key value pairs of element symbols
            and their log epsilon abundances
        '''

        elem_num_dict = dict(sorted({atomic_sym_to_num(
            elem): val for elem, val in abund_dict.items()}.items()))
        with open(self.file_loc, 'r') as file:
            lines = file.readlines()

        newfilename = 'new_atmos_file.tmp'
        with open(newfilename, 'w') as newfile:
            for line in lines:
                if 'NATOMS' in line:
                    nline = line.replace('NATOMS     0',
                                         f'NATOMS    {len(elem_num_dict):2d}')
                    print(nline, file=newfile, end='')
                    for elem, val in elem_num_dict.items():
                        print(f'   {elem}  {val:.2f}', file=newfile)
                else:
                    print(line, file=newfile, end='')

        os.replace(newfilename, self.file_loc)
