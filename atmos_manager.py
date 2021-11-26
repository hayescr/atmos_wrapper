import os
from abund_utils import atomic_sym_to_num


class AtmosManager:

    wrapper_path = os.path.dirname(os.path.realpath(__file__))

    def __init__(self, teff, logg, feh, vmicro=None, cfe=0., alphafe=0.):
        # Need to decide if we are going to have c and alpha input as flags
        # Or if we want to just input values and just force add set C and alpha
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

    def moog_abund(self, abund_dict):
        # To be implemented
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
                        print(f'   {elem}  {val}', file=newfile)
                else:
                    print(line, file=newfile, end='')

        os.replace(newfilename, self.file_loc)
