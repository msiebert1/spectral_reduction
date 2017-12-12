
def reduce(imglist, _listsens, _listarc, _ext_trace, _dispersionline, _cosmic, _interactive):
    # print "LOGX:: Entering `efoscfastredu` method/function in %(__file__)s"
    # % globals()

    import string
    import os
    import re
    import sys
    os.environ["PYRAF_BETA_STATUS"] = "1"
    try:      from astropy.io import fits as pyfits
    except:   import   pyfits
    # from ntt.util import readhdr, readkey3
    # import ntt
    import numpy as np
    import util
    import instruments
    import combine_sides as cs
    import sys
    sys.path.append("../cosmics.py_0.4/")
    import cosmics
    dv = util.dvex()
    scal = np.pi / 180.
    if not _interactive:
        _interactive = False
        _inter = 'NO'
    else:
        _inter = 'YES'
    from pyraf import iraf

    hdr0 = util.readhdr(imglist[0])
    if util.readkey3(hdr0, 'VERSION') == 'kastb':
        inst = instruments.kast_blue
    elif util.readkey3(hdr0, 'VERSION') == 'kastr':
        inst = instruments.kast_red
    else:
        print util.readkey3(hdr0, 'VERSION') + 'not in database'
        sys.exit()

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['ccdproc', 'imcopy', 'specred.apall', 'longslit.identify', 'longslit.reidentify', 'specred.standard',
                'longslit.fitcoords', 'onedspec.wspectext']
    for t in toforget:
        iraf.unlearn(t)
    iraf.ccdred.verbose = 'no'  # not print steps
    iraf.specred.verbose = 'no'  # not print steps
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.ccdtype = ''

    iraf.longslit.mode = 'h'
    iraf.specred.mode = 'h'
    iraf.noao.mode = 'h'
    iraf.ccdred.instrument = "ccddb$kpno/camera.dat"
    # iraf.set(direc=ntt.__path__[0] + '/')

    asci_files = []
    newlist = [[],[]]

    if _cosmic:
        for i in imglist:
            print i
            if 'b' in i:
                inst = instruments.kast_blue
            elif 'r' in i:
                inst = instruments.kast_red
            array, header = cosmics.fromfits(i)
            c = cosmics.cosmicsimage(array, gain=inst.get('gain'), readnoise=inst.get('read_noise'), sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)
            c.run(maxiter = 4)
            cosmics.tofits('c' + i, c.cleanarray, header)
        for j in range(len(imglist)):
            imglist[j] = 'c' + imglist[j]

    print imglist

    for img in imglist:
        if 'b' in img:
            newlist[0].append(img)
        elif 'r' in img:
            newlist[1].append(img)

    if len(newlist[1]) < 1:
        newlist = newlist[:-1]
    # for img in imglist:
    print newlist
    for imgs in newlist:
        hdr = util.readhdr(imgs[0])
        # hdr = util.readhdr(img)
        # _tech = util.readkey3(hdr, 'tech')

        # print hdr
        if util.readkey3(hdr, 'VERSION') == 'kastb':
            inst = instruments.kast_blue
        elif util.readkey3(hdr, 'VERSION') == 'kastr':
            inst = instruments.kast_red
        else:
            print util.readkey3(hdr, 'VERSION') + 'not in database'
            sys.exit()

        iraf.specred.dispaxi = inst.get('dispaxis')
        iraf.longslit.dispaxi = inst.get('dispaxis')

        _gain = inst.get('gain')
        _ron = inst.get('read_noise')
        iraf.specred.apall.readnoi = _ron
        iraf.specred.apall.gain = _gain

        print inst.get('name')
        print imgs

        # if _tech != 'SPECTRUM':
        #     sys.exit('error: ' + str(img) + ' is not a spectrum ')
        # print '\n####  image name = ' + img + '\n'
        # _grism0 = inst.get('grism')
        # _filter0 = inst.get('filter')
        # _slit0 = inst.get('slit')
        _object0 = util.readkey3(hdr, 'OBJECT')
        _date0 = util.readkey3(hdr, 'DATE-OBS')

        # setup = (_grism0, _filter0, _slit0)
        #TODO: change these for different setups
        _biassec0 = inst.get('biassec')
        _trimsec0 = inst.get('trimsec')

        # if _grism0 == 'Gr16':
        #     _trimsec0 = '[100:950,1:950]'
        # elif _grism0 == 'Gr13':
        #     if _filter0 == 'Free':
        #         _trimsec0 = '[100:950,1:1015]'
        #     elif _filter0 == 'GG495':
        #         _trimsec0 = '[100:950,208:1015]'
        #     elif _filter0 == 'OG530':
        #         _trimsec0 = '[100:950,300:1015]'
        # elif _grism0 == 'Gr11':
        #     _trimsec0 = '[100:950,5:1015]'
        # else:
        #     _trimsec0 = '[100:950,5:1015]'

        _object0 = re.sub(' ', '', _object0)
        _object0 = re.sub('/', '_', _object0)
        nameout0 = str(_object0) + '_' + inst.get('name') + '_' + str(_date0)
        # for _set in setup:
        #     nameout0 = nameout0 + '_' + _set
        # nameout0 = util.name_duplicate(img, nameout0, '')
        nameout0 = util.name_duplicate(imgs[0], nameout0, '')
        timg = nameout0

        if len(imgs) > 1:
            img_str = ''
            for i in imgs:
                img_str = img_str + i + ','
            iraf.imcombine(img_str, output=timg)
        else:
            img = imgs[0]
            if os.path.isfile(timg):
                os.system('rm -rf ' + timg)
            iraf.imcopy(img, output=timg)
        print img
        
        
        #TODO: Update with masterbias/master flat
        iraf.ccdproc(timg, output='', overscan='yes', trim='yes', zerocor="no", flatcor="no", readaxi='line',
                     trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)
        
        img = timg

        if _listarc:
            arcfile = util.searcharc(img, _listarc)[0]
        else:
            arcfile = ''
        if not arcfile:
            arcname = raw_input('Arc Filename: ')
            arcfile = arcname
            util.delete('t' + arcfile)
            iraf.ccdproc(arcfile, output= 't' + arcfile, overscan='yes', trim='yes', zerocor="no", flatcor="no",
                         readaxi='line', trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)
            arcfile = 't' + arcfile
            # raise TypeError
        else:
            #TODO: Update with masterbias/master flat
            iraf.ccdproc(arcfile, output='t' + arcfile, overscan='yes', trim='yes', zerocor="no", flatcor="no",
                         readaxi='line', trimsec=str(_trimsec0), biassec=str(_biassec0), Stdout=1)
            arcfile = 't' + arcfile

        # if _cosmic:
        #     # print cosmic rays rejection
        #     ntt.cosmics.lacos(img, output='', gain=_gain, readn=_ron, xorder=9, yorder=9, sigclip=4.5, sigfrac=0.5,
        #                       objlim=1, verbose=True, interactive=False)
        #     print '\n### cosmic rays rejections ........ done '

        if not arcfile:
            print '\n### warning no arcfile \n exit '
        else:
            arcref = inst.get('QL_arc')
            if arcfile[0] == '/':
                os.system('cp ' + arcfile + ' ' +
                          string.split(arcfile, '/')[-1])
                arcfile = string.split(arcfile, '/')[-1]
            # arcref = string.split(arcref, '/')[-1]
            if arcref:
                os.system('cp ' + arcref + ' .')
                # arcref = string.split(arcref, '/')[-1]
                if not os.path.isdir('database/'):
                    os.mkdir('database/')
                # if os.path.isfile(util.searcharc(img, '')[1] + '/database/id' + re.sub('.fits', '', arcref)):
                #     os.system('cp ' + util.searcharc(img, '')[1] + '/database/id' + re.sub('.fits', '',
                #                                                                                arcref) + ' database/')
                #TODO: Generalize this path
                # os.system('cp ' + arcref + ' /home/msiebert/Documents/UCSC/Research/spectral_reduction/database/id' + re.sub('.fits', '', arcref))
                
                print arcref, arcfile
                # if inst.get('rotate'):
                #     arcref = arcref + '[0]'
                #     arcfile = arcfile + '[0]'
                    # iraf.imtranspose(input=arcref + '[-*,*]', output=re.sub('.fits', '', arcref)+'rot')
                    # iraf.imtranspose(input=arcfile + '[-*,*]', output=re.sub('.fits', '', arcfile)+'rot')
                    # iraf.imtranspose(input=img + '[-*,*]', output=re.sub('.fits', '', img)+'rot')
                    # arcref = re.sub('.fits', '', arcref) + 'rot'
                # print inst.get('ident_order')
                # iraf.longslit.identify(images=arcref, section=inst.get('section'),
                #                        coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', nsum=10, fwidth=7,
                #                        order=3, mode='h')
                # iraf.longslit.identify(images=arcref, section=inst.get('section'),
                #                        coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', nsum=10, fwidth=7,
                #                        order=4, function='legendre', mode='h')
                # raise TypeError
                # arcref = re.sub('.fits', '', arcref)
                # arcfile = re.sub('.fits', '', arcfile)
                _inter='NO'
                iraf.longslit.reidentify(referenc=arcref, images=arcfile, interac=_inter, section=inst.get('section'),
                                         coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', overrid='yes', step=0,
                                         newaps='no', nsum=5, nlost=2, mode='h', verbose='no')
                # raise TypeError
            else:
                iraf.longslit.identify(images=arcfile, section=inst.get('section'),
                                       coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', nsum=10, fwidth=7,
                                       order=3, mode='h')

            iraf.longslit.reident(referenc=arcfile, images=arcfile, interac='NO', section=inst.get('section'),
                                  coordli='direc$standard/ident/Lines_HgCdHeNeAr600.dat', overrid='yes', step=10,
                                  newaps='yes', nsum=5, nlost=2, mode='h', verbose='no')
            # qqq = iraf.longslit.fitcoords(images=re.sub('.fits', '', arcfile), fitname= re.sub('.fits', '', arcfile),
            #                               interac='no', combine='yes', databas='database',
            #                               function='legendre', yorder=4, logfile='logfile', plotfil='', mode='h')
            qqq = iraf.longslit.fitcoords(images=re.sub('.fits', '', arcfile), fitname= re.sub('.fits', '', arcfile),
                                          interac='no', combine='yes', databas='database',
                                          function='legendre', yorder=4, logfile='logfile', plotfil='', mode='h')

            print img, arcfile, arcref
            iraf.specred.transform(input=img, output=img, minput='', fitnames=re.sub('.fits', '', arcfile),
                                   databas='database',
                                   x1='INDEF', x2='INDEF', y1='INDEF', y2='INDEF', flux='yes', mode='h',
                                   logfile='logfile')

            # raise TypeError
            # ######################  check wavelength calibration ############
            # _skyfile = ntt.__path__[
            #     0] + '/standard/ident/sky_' + setup[0] + '_' + setup[1] + '.fits'
            # shift = ntt.efoscspec2Ddef.skyfrom2d(img, _skyfile)
            # print '\n###     check in wavelengh performed ...... spectrum shifted of  ' + str(shift) + ' Angstrom \n'
            # zro = pyfits.open(img)[0].header.get('CRVAL2')
            # ntt.util.updateheader(img, 0, {'CRVAL2': [zro + int(shift), '']})
            # std, rastd, decstd, magstd = ntt.util.readstandard(
            #     'standard_efosc_mab.txt')
            # hdrt = readhdr(img)
            # _ra = readkey3(hdrt, 'RA')
            # _dec = readkey3(hdrt, 'DEC')
            # _object = readkey3(hdrt, 'object')
            # dd = np.arccos(np.sin(_dec * scal) * np.sin(decstd * scal) + np.cos(_dec * scal) *
            #                np.cos(decstd * scal) * np.cos((_ra - rastd) * scal)) * ((180 / np.pi) * 3600)
            # if min(dd) < 100:
            #     _type = 'stdsens'
            #     ntt.util.updateheader(
            #         img, 0, {'stdname': [std[np.argmin(dd)], '']})
            #     ntt.util.updateheader(
            #         img, 0, {'magstd': [float(magstd[np.argmin(dd)]), '']})
            # else:
            #     _type = 'obj'
            print '\n###      EXTRACTION USING IRAF TASK APALL \n'
            result = []
            _type ='obj'
            if _type == 'obj':
                # _interactive = 'NO'
                imgex = util.extractspectrum(
                    img, dv, inst, _ext_trace, _dispersionline, _interactive, _type)
                # raise TypeError
                #TEST STANDARD STAR
                # standardfile = imgex
                # _airmass = readkey3(hdr, 'AIRMASS')
                # _exptime = readkey3(hdr, 'EXPTIME')
                # _caldir = 'direc$standard/MAB/'
                # _extinctdir = 'direc$standard/extinction/'
                # refstar =  'm' + re.sub('.dat', '', 'feige34')

                # iraf.specred.standard(input=standardfile, output='std_feige34.fits',
                #                   caldir=_caldir, star_nam=refstar, airmass=_airmass,
                #                   exptime=_exptime, interac=_interactive)
                # iraf.specred.sensfunc(standard='std_feige34.fits', sensitiv='sens_feige34.fits',
                #                   ignorea='yes', functio='spline3', order=16,
                #                   interac=_interactive)
                # util.updateheader(
                #     imgex, 0, {'FILETYPE': [22107, 'extracted 1D spectrum ']})
                # util.updateheader(imgex, 0, {
                #     'PRODCATG': ['SCIENCE.' + readkey3(readhdr(imgex), 'tech').upper(), 'Data product category']})
                # util.updateheader(imgex, 0, {'TRACE1': [img, '']})
                # result.append(imgex)
                # raise TypeError
                sensfile = inst.get('QL_sens')
                if sensfile:
                    print sensfile
                    imgf = re.sub('.fits', '_f.fits', img)
                    _extinctdir = 'direc$standard/extinction/'
                    _extinction = 'extinction_lasilla.dat'
                    _observatory = 'lasilla'
                    _exptime = util.readkey3(hdr, 'EXPTIME')
                    _airmass = util.readkey3(hdr, 'AIRMASS')
                    util.delete(imgf)
                    # iraf.specred.calibrate(input=imgex, output=imgf, sensiti=sensfile, extinct='yes',
                    #                        flux='yes', ignorea='yes', extinction=_extinctdir + _extinction,
                    #                        observatory=_observatory, airmass=_airmass, exptime=_exptime,
                    #                        fnu='no')
                    iraf.specred.calibrate(input=imgex, output=imgf, sensiti=sensfile, extinct='no',
                                           flux='yes', ignorea='yes', airmass=_airmass, exptime=_exptime,
                                           fnu='no')
                    # hedvec = {'SENSFUN': [string.split(sensfile, '/')[-1], 'sensitivity function'],
                    #           'FILETYPE': [22208, '1D wavelength and flux calibrated spectrum '],
                    #           'SNR': [util.StoN2(imgf, False), 'Average S/N ratio'],
                    #           'BUNIT': ['erg/cm2/s/Angstrom', 'Flux Calibration Units'], 'TRACE1': [imgex, '']}
                    # util.updateheader(imgf, 0, hedvec)
                    imgout = imgf
                    # imgd = efoscspec1Ddef.fluxcalib2d(img, sensfile)
                    # util.updateheader(
                    #     imgd, 0, {'FILETYPE': [22209, '2D wavelength and flux calibrated spectrum ']})
                    # util.updateheader(imgd, 0, {'TRACE1': [img, '']})
                    imgasci = re.sub('.fits', '.asci', imgout)
                    errasci = re.sub('.fits', '_err.asci', imgout)
                    util.delete(imgasci)
                    iraf.onedspec.wspectext(
                        imgout + '[*,1,1]', imgasci, header='no')
                    iraf.onedspec.wspectext(
                        imgout + '[*,1,4]', errasci, header='no')
                    spec = np.transpose(np.genfromtxt(imgasci))
                    err = np.transpose(np.genfromtxt(errasci))
                    util.delete(errasci)
                    final = np.transpose([spec[0], spec[1], err[1]])
                    np.savetxt(imgasci, final)

                    # result = result + [imgout, imgd, imgasci]
                    result = result + [imgout, imgasci]
            else:
                imgex = util.extractspectrum(
                    img, dv, _ext_trace, _dispersionline, _interactive, 'std')
                imgout = efoscspec1Ddef.sensfunction(
                    imgex, 'spline3', 6, _inter)
                result = result + [imgout]

    # for img in result:
    #     if img[-5:] == '.fits':
    #         ntt.util.phase3header(img)  # phase 3 definitions
    #         ntt.util.airmass(img)  # phase 3 definitions
    #         ntt.util.updateheader(
    #             img, 0, {'quality': ['Rapid', 'Final or Rapid reduction']})

        result = result + [imgex] + [timg]
        print result
        asci_files.append(imgasci)
        if not os.path.isdir(_object0 + '/'):
            os.mkdir(_object0 + '/')
            for img in result:
                os.system('mv ' + img + ' ' + _object0 + '/')
        else:
            for img in result:
                os.system('mv ' + img + ' ' + _object0 + '/')


    if len(asci_files) > 1:
        final = cs.combine_blue_red(asci_files[0], asci_files[1], _object0)
    return result