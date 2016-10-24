__author__ = 'ltesti - 8 Jul 2015'

import numpy as np
from star_data import StarData
import scipy.constants as sc
import scipy.stats as ss
import os
import matplotlib.pyplot as plt
from statsmodels.nonparametric import smoothers_lowess as sml


class UVData(object):
    def __init__(self, star, file_tuple, doread=True):
        if doread:
            self.infile = file_tuple[0]
            self.wl = file_tuple[1]
            self.__read_uv_from_ascii_file()
        self.__set_star_infile(star)
        self.__calc_uvd_amp_dep(in_klambda=True)

    def __calc_uvd_amp_dep(self, in_klambda=True):
        self.uvd = self.calc_uvd(self.u, self.v, in_klambda=in_klambda)
        self.amp = self.calc_amp(self.re, self.im)
        self.deproject(self.myStar.pa, self.myStar.inc)

    def __set_star_infile(self, sname):
        self.myStar = StarData(sname)

    def __str__(self):
        rep = "UVData for: " + self.myStar.name + ", read from ascii file: " + self.infile + ", wavelength: " + str(
            self.wl)

        return rep

    # reads from standard uvfile from disk (u,v,re,im,w)
    def __read_uv_from_ascii_file(self):
        openfile = open(self.infile, 'r')
        lines = openfile.readlines()
        openfile.close()

        size = len(lines)
        self.u = np.zeros(size, dtype='float64')
        self.v = np.zeros(size, dtype='float64')
        self.re = np.zeros(size, dtype='float64')
        self.im = np.zeros(size, dtype='float64')
        self.we = np.zeros(size, dtype='float64')
        for j, line in enumerate(lines):
            self.u[j], self.v[j], self.re[j], self.im[j], self.we[j] = line.split()[0:5]

            # self.wl = intuple[1]

    # Compute/set uvdistance, wle in the same units as (u,v)
    def calc_uvd(self, u, v, in_klambda=True):
        factor = 1.
        if in_klambda:
            factor = 1.e3
        return np.sqrt(u * u + v * v) / self.wl / factor

    # Compute/set amplitude
    @staticmethod
    def calc_amp(re, im):
        return np.sqrt(re * re + im * im)

    # Compute average and dispersion of re and im
    # note: here re, im, we are arbitrary vecors, I did
    #       it this way to allow easy means on subsets
    #       another (better?) option would be to pass the
    #       indices of the re, im, we that we want to use
    @staticmethod
    def calc_avg(re, im, we):
        ws = we.sum()
        if ws > 0.:
            rem = (re * we).sum() / ws
            imm = (im * we).sum() / ws
        else:
            rem = 0.
            imm = 0.
        erem = np.sqrt(((re - rem) * (re - rem)).sum()) / (len(re) - 1.)
        eimm = np.sqrt(((im - imm) * (im - imm)).sum()) / (len(im) - 1.)

        return rem, erem, imm, eimm

    # Compute binned re, im and amp
    # returns the binned arrays and std deviation in each bin
    def calc_binned(self, bin, dep=True):
        if dep:
            uvd = self.uvd_dep
        else:
            uvd = self.uvd
        db2 = (bin[1] - bin[0]) / 2.
        reb = np.zeros(len(bin))
        imb = np.zeros(len(bin))
        amb = np.zeros(len(bin))
        ereb = np.zeros(len(bin))
        eimb = np.zeros(len(bin))
        eamb = np.zeros(len(bin))
        for i in xrange(len(bin)):
            idx = np.where((uvd > bin[i] - db2) & (uvd <= bin[i] + db2))
            reb[i], ereb[i], imb[i], eimb[i] = self.calc_avg(self.re[idx], self.im[idx], self.we[idx])
            amb[i] = np.sqrt(reb[i] * reb[i] + imb[i] * imb[i])
            eamb[i] = np.sqrt(((amb[i]-amb[idx])**2).sum())/(len(idx)-1.)
            eamb[i] = np.sqrt(ereb[i] * ereb[i] + eimb[i] * eimb[i])

        uv_binned = dict(reb=reb, ereb=ereb, imb=imb, eimb=eimb, amb=amb, eamb=eamb)
        return uv_binned

    # Deproject visibilities for a position angle pa and inclination inc
    def deproject(self, pa, inc):

        if inc == 0.0:
            self.u_dep = np.copy(self.u)
            self.v_dep = np.copy(self.v)
            self.uvd_dep = np.copy(self.uvd)
        else:
            ur = self.u * np.cos(pa) - self.v * np.sin(pa)
            vr = self.u * np.sin(pa) + self.v * np.cos(pa)

            ur *= np.cos(inc)

            self.u_dep = ur * np.cos(-pa) - vr * np.sin(-pa)
            self.v_dep = ur * np.sin(-pa) + vr * np.cos(-pa)

            self.uvd_dep = self.calc_uvd(self.u_dep, self.v_dep, in_klambda=True)

    def uvplot(self, yp="Real", uv_dep=False):
        plot_names = dict(real='real', re='real',
                          imag='imag', imaginary='imag', im='imag',
                          amp='amp', amplitude='amp')
        uvp = self.uvd
        if uv_dep:
            uvp = self.uvd_dep
        if yp.lower() in plot_names:
            doplot = plot_names[yp.lower()]
            if doplot == "real":
                yp = self.re
            if doplot == "imag":
                yp = self.im
            if doplot == "amp":
                yp = self.amp
        else:
            printf("Error! Variable {0} not recognized!".format(yp))
            printf("Plotting amplitude...")
            yp = self.amp
        plt.plot(uvp, yp, ".", color="blue")

    @staticmethod
    def __fg(x, s):
        ff = np.exp(-x * x / (2. * s * s))
        return ff

    @staticmethod
    def __get_fit(x, y, my_is_sorted=True, my_delta=4., my_frac=0.05, my_it=3, my_return_sorted=False):
        y_f = sml.lowess(y, x, is_sorted=my_is_sorted, delta=my_delta, frac=my_frac, it=my_it, return_sorted=my_return_sorted)
        # if plotfit:
        #     plt.plot(y_f[:, 0], y_f[:, 1], color='green', linestyle='dashed', linewidth=2)
        return y_f

    def __do_plot_weights(self, xp, re_f, im_f, re_r, re_std, im_r, im_std, plsig=1.6):
        #
        # Plot real and imaginary part
        # get uvd minmax for plot
        duvd = max(self.uvd) - min(self.uvd)
        uvp_min = min(self.uvd) - 0.025*duvd
        uvp_max = max(self.uvd) + 0.025*duvd
        #
        # Top panel is Real part uvplot
        plt.subplot(4, 1, 1)
        self.uvplot(yp='real', uv_dep=False)
        plt.plot(xp,re_f,color='green',linestyle='dashed',linewidth=2)
        plt.plot([0., 1.e9], [0., 0.], color='green', linestyle='dotted')
        plt.xlim(uvp_min, uvp_max)
        #
        # Middle panel is Imaginary part uvplot
        plt.subplot(4, 1, 2)
        self.uvplot(yp='imag', uv_dep=False)
        plt.plot(xp, im_f, color='green', linestyle='dashed', linewidth=2)
        plt.plot([0., 1.e9], [0., 0.], color='green', linestyle='dotted')
        plt.xlim(uvp_min, uvp_max)
        #
        # Bottom left is Real part histogram
        plt.subplot(2, 2, 3)
        self.__subplot_hist_gauss(re_r, re_std, pl_sig=plsig, mycolor='red')
        #
        # Bottom left is Imaginary part histogram
        plt.subplot(2, 2, 4)
        self.__subplot_hist_gauss(im_r, im_std, pl_sig=plsig, mycolor='blue')

    def __subplot_hist_gauss(self, v, std, pl_sig=1.6, mycolor='red'):
            dxp = pl_sig*std
            nn, xb, patches = plt.hist(v, 100, color=mycolor)
            plt.xlim(-dxp,dxp)
            xg = np.linspace(-dxp,dxp, 100)
            yy = self.__fg(xg, std)
            ff = (nn.sum() * (xb[1] - xb[0])) / (yy.sum() * (xg[1] - xg[0]))
            plt.plot(xg, ff * yy, linestyle="dashed", color='green', linewidth=2)

    #
    # Calcola la standard deviation usando un sigma clip
    #    sulla prima iterazione
    @staticmethod
    def __get_std(v, sigma_clip=2.0):
        std_0 = ss.tstd(v)
        dv = sigma_clip * std_0
        return ss.tstd(v, limits=(-dv, dv))

    def get_weight(self, doplot=True):
        #
        # Sort vectors
        re_sort = np.copy(self.re[np.argsort(self.uvd)])
        im_sort = np.copy(self.im[np.argsort(self.uvd)])
        uvd_sort = np.copy(self.uvd[np.argsort(self.uvd)])
        #
        # Fit and subtract real and imaginary part
        delta_uvd = 0.01*(uvd_sort.max()-uvd_sort.min())
        re_f = self.__get_fit(uvd_sort, re_sort, my_delta=delta_uvd)
        re_r = re_sort - re_f
        im_f = self.__get_fit(uvd_sort, im_sort, my_delta=delta_uvd)
        im_r = im_sort - im_f
        #
        #
        re_std = self.__get_std(re_r, sigma_clip=1.55)
        im_std = self.__get_std(im_r, sigma_clip=1.55)
        #
        # Plot proudly, if requested...
        if doplot:
            self.__do_plot_weights(uvd_sort, re_f, im_f, re_r, re_std, im_r, im_std, plsig=3.0)
        #
        # Calcolo dei fattori per i pesi
        wm = np.sum(self.we) / len(self.we)
        rat_re = 1. / (wm * re_std ** 2.)
        rat_im = 1. / (wm * im_std ** 2.)
        print("Ratios (1./std^2)/we: Real={0}  Imag={1}  Average={2}".format(rat_re, rat_im, (rat_re + rat_im) / 2.))

        return rat_re, rat_im

    def write_uv_to_ascii(self,outfile):
        #
        # Open the output file (Format will be: u,v,Re,Im,We,spw)
        f_out = open(outfile, "w")

        for i in range(len(self.u)):
            f_out.write("%+.15e  %+.15e  %+.10e  %+.10e  %+.10e\n" % (self.u[i],self.v[i],self.re[i],self.im[i],self.we[i]))
            #print >> f_out, "%+.15e  %+.15e  %+.10e  %+.10e  %+.10e  %d" % \
            #         (self.u[i],self.v[i],self.re[i],self.im[i],self.we[i],self.spw[i])

        # Close output file
        f_out.close()



class UVDataMS(UVData):
    # Modifies the class UVData to read the data from a MS file
    # required to be in CASA and have access to the tb tool

    def __init__(self, star, file_tuple):
        self.infile = file_tuple[0]
        self.tb = file_tuple[1]
        self.dualpol = True
        if len(file_tuple) > 2:
            self.dualpol = file_tuple[2]
        self.__read_uv_from_ms_file()
        super(UVDataMS, self).__init__(star, file_tuple, doread=False)

    def __str__(self):
        rep = "UVData for: " + self.myStar.name + ", read from ms file: " + self.infile + ", wavelength: " + str(
            self.wl)

        return rep

    #
    # reads uvdata from MS file
    #
    def __read_uv_from_ms_file(self):

        #
        # Read columns from the ms file
        data, uvw, weight, spw, spec_tab = self.__read_ms_table()

        #
        # Copy over the relevant column data into the appropriate np arrays
        # Takes care of dualpol data: note the two pols are averaged (weighted)
        self.u = np.array(uvw[0, :])
        self.v = np.array(uvw[1, :])
        self.w = np.array(uvw[2, :])
        re_xx = np.array(data[0, 0, :].real)
        im_xx = np.array(data[0, 0, :].imag)
        wei_xx = np.array(weight[0, :])
        if self.dualpol:
            re_yy = data[1, 0, :].real
            im_yy = data[1, 0, :].imag
            wei_yy = weight[1, :]
        else:
            re_yy = data[0, 0, :].real
            im_yy = data[0, 0, :].imag
            wei_yy = weight[0, :]
        # Here we actually write the real, imaginary and weights into the arrays
        # Note that we have a safeguard for flagged data with we=0
        self.we = (wei_xx + wei_yy)/2.
        self.re = np.copy(self.we)
        self.im = np.copy(self.we)
        nw_good = np.where(self.we != 0.0)
        self.re[nw_good] = (re_xx[nw_good] * wei_xx[nw_good] + re_yy[nw_good] * wei_yy[nw_good]) / (
            wei_xx[nw_good] + wei_yy[nw_good])
        self.im[nw_good] = (im_xx[nw_good] * wei_xx[nw_good] + im_yy[nw_good] * wei_yy[nw_good]) / (
            wei_xx[nw_good] + wei_yy[nw_good])

        # From the SPECTRAL_WINDOW table extract channel frequencies from CHAN_FREQ and
        # compute average freq, then convert to wavelength in m
        self.tb.open(spec_tab)
        self.wl = sc.c / self.tb.getcol('CHAN_FREQ').mean()
        self.tb.close()

    #
    # This is the function that actually reads the columns from the MS file on disk
    def __read_ms_table(self):
        # Open ms table
        self.tb.open(self.infile)

        # Check if CORRECTED_DATA column is present:
        all_columns = self.tb.colnames()
        self.correct = False
        if 'CORRECTED_DATA' in all_columns:
            data = self.tb.getcol("CORRECTED_DATA")
            self.correct = True
        else:
            data = self.tb.getcol("DATA")

        # Get column UVW, WEIGHT and DATA_DESC_ID (spw info!)
        uvw = self.tb.getcol("UVW")
        weight = self.tb.getcol("WEIGHT")
        spw = self.tb.getcol("DATA_DESC_ID")
        spec_tab = self.tb.getkeyword('SPECTRAL_WINDOW')[7:]

        # Close ms file
        self.tb.close()

        return data, uvw, weight, spw, spec_tab

    #
    # write uv-data to MS file using as template the
    # input file for the class definition
    #
    def write_uv_to_ms_file(self, new_ms):
        #
        # Create new table copying the reference one
        syscommand = 'rm -rf ' + new_ms
        os.system(syscommand)
        syscommand = 'cp -r ' + self.infile + ' ' + new_ms
        os.system(syscommand)

        # Open the new MS, read the current data from the reference table
        data, uvw, weight, spw, spec_tab = self.__read_ms_table()

        #
        # Prepare columns with new data to be written in the table
        data[0, 0, :] = self.re + 1j*self.im
        weight[0, :] = self.we
        if self.dualpol:
            data[1, 0, :] = self.re + 1j*self.im
            weight[1, :] = self.we

        #
        # Write out the new data in the (copied) table columns
        self.tb.open(new_ms, nomodify=False)
        self.tb.putcol("DATA", data)
        self.tb.putcol("WEIGHT", weight)
        if self.correct is True:
            self.tb.putcol("CORRECTED_DATA", data)
        self.tb.flush()
        # Close table tools
        self.tb.close()
