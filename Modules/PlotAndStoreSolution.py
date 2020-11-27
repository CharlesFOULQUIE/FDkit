#!/usr/bin/env python

import time, glob, shutil, os
import numpy as np

def merge_zip_archives(individual_archives, archive_name):
    """
    Merge individual zip archives made with numpy.savez into
    one archive with name archive_name.
    The individual archives can be given as a list of names
    or as a Unix wild chard filename expression for glob.glob.
    The result of this function is that all the individual
    archives are deleted and the new single archive made.
    """
    import zipfile
    archive = zipfile.ZipFile(
        archive_name, 'w', zipfile.ZIP_DEFLATED,
        allowZip64=True)
    if isinstance(individual_archives, (list,tuple)):
        filenames = individual_archives
    elif isinstance(individual_archives, str):
        filenames = glob.glob(individual_archives)

    # Open each archive and write to the common archive
    for filename in filenames:
        f = zipfile.ZipFile(filename,  'r',
                            zipfile.ZIP_DEFLATED)
        for name in f.namelist():
            data = f.open(name, 'r')
            # Save under name without .npy
            archive.writestr(name[:-4], data.read())
        f.close()   
    archive.close()
    


class PlotAndStoreSolution:
    """
    Class for the user_action function in solver.
    Visualizes the solution only.
    """
    def __init__(
        self,
        casename='tmp',    # Prefix in filenames
        umin=-1, umax=1,   # Fixed range of y axis
        pause_between_frames=None,  # Movie speed
        screen_movie=True, # Show movie on screen?
        title='',          # Extra message in title
        skip=1,            # Skip Plot and stor solution
        filename=None,     # Name of file with solutions 
        init=None):
        self.res = 'res'
        self.casename = casename
        self.yaxis = [umin, umax]
        self.pause = pause_between_frames
        import matplotlib.pyplot as plt
        self.plt = plt
        self.screen_movie = screen_movie
        self.title = title
        self.skip = skip
        self.filename = filename
        self.init = init
        if filename is not None:
            # Store time points when u is written to file
            self.t = []
            filenames = glob.glob(os.path.join(self.res,'.' + self.filename + '*.dat.npz'))
            for filename in filenames:
                os.remove(filename)

        # Clean up old movie frames
        for filename in glob.glob(os.path.join(self.res,'frame_*.png')):
            os.remove(filename)

    def __call__(self, u, x, t, n):
        """
        Callback function user_action, call by solver:
        Store solution, plot on screen and save to file.
        """
        # Save archive
        if n % self.skip != 0:
            return
            
        # Store solution
        if self.filename is not None:
            name = 'u%04d' % n  # array name
            kwargs = {name: u}
            fname = '.' + self.filename + '_' + name + '.dat'
            np.savez(os.path.join(self.res,fname), **kwargs)
            self.t.append(t[n])  # store corresponding time value
            if n == 0:           # save x once
                np.savez(os.path.join(self.res,'.' + self.filename + '_x.dat'), x=x)
                
        # Plote solution     
        if self.screen_movie is not False :
    
            umin, umax = self.yaxis
            title = 'Nx=%d' % (x.size-1)
            if self.title:
                title = self.title + ' ' + title

            # native matplotlib animation
            if n == 0:
                self.plt.ion()
                self.lines = self.plt.plot(x, self.init(x), 'k--')
                self.lines = self.plt.plot(x, u, 'b-')
                self.plt.axis([x[0], x[-1], umin, umax])
                self.plt.xlabel('x')
                self.plt.ylabel('u')
                self.plt.title(title)
                self.plt.legend(['t=%.3f' % t[n]])
                self.plt.pause(0.0001) 
            else:
                # Update new solution
                self.lines[0].set_ydata(u)
                self.plt.legend(['t=%.3f' % t[n]])
                self.plt.draw()
                self.plt.pause(0.0001) 
            
            # pause
            if t[n] == 0:
                time.sleep(2)  # let initial condition stay 2 s
            else:
                if self.pause is None:
                    pause = 0.0001 if u.size < 100 else 0
                time.sleep(pause)
    
            self.plt.savefig(os.path.join(self.res,'frame_%04d.png' % (n/self.skip)))

            if n == (len(t) - 1):   # finished with this run, close plot
                self.plt.close()
        else:
            print(n)

    def make_movie_file(self):
        """
        Create subdirectory based on casename, move all plot
        frame files to this directory, and generate
        an index.html for viewing the movie in a browser
        (as a sequence of PNG files).
        """
        os.chdir(self.res)
        # Make HTML movie in a subdirectory
        directory = self.casename

        if os.path.isdir(directory):
            shutil.rmtree(directory)   # rm -rf directory
        os.mkdir(directory)            # mkdir directory
        # mv frame_*.png directory
        for filename in glob.glob('frame_*.png'):
            os.rename(filename, os.path.join(directory, filename))
        os.chdir(directory)        # cd directory

        fps = 12 # frames per second

        # Make other movie formats: Flash, Webm, Ogg, MP4
        codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                         libtheora='ogg')
        filespec = 'frame_%04d.png'
        movie_program = 'ffmpeg'  # or 'avconv'
        for codec in codec2ext:
            ext = codec2ext[codec]
            cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
                  '-vcodec %(codec)s movie.%(ext)s' % vars()
            os.system(cmd)

        os.chdir(os.pardir)  # move back to parent directory
        os.chdir(os.pardir)

    def get_res_directory(self):
        return self.res
        
    def close_file(self, hashed_input):
        """
        Merge all files from savez calls into one archive.
        hashed_input is a string reflecting input data
        for this simulation (made by solver).
        """
        if self.filename is not None:
            # Save all the time points where solutions are saved
            np.savez(os.path.join(self.res,'.' + self.filename + '_t.dat'),
                     t=np.array(self.t, dtype=float))

            # Merge all savez files to one zip archive
            archive_name = os.path.join(self.res,'.' + hashed_input + '_archive.npz')
            filenames = glob.glob(os.path.join(self.res,'.' + self.filename + '*.dat.npz'))
            merge_zip_archives(filenames, archive_name)
            print('Archive name:', archive_name)
            
            # Suppress tmp archive
            filenames = glob.glob(os.path.join(self.res,'.' + self.filename + '*.dat.npz'))
            for filename in filenames:
                os.remove(filename)

class PlotMediumAndSolution(PlotAndStoreSolution):
    def __init__(self, medium, **kwargs):
        """Mark medium in plot: medium=[x_L, x_R]."""
        self.medium = medium
        PlotAndStoreSolution.__init__(self, **kwargs)

    def __call__(self, u, x, t, n):

        # Save archive
        if n % self.skip != 0:
            return
            
        # Store solution
        if self.filename is not None:
            name = 'u%04d' % n  # array name
            kwargs = {name: u}
            fname = '.' + self.filename + '_' + name + '.dat'
            np.savez(os.path.join(self.res,fname), **kwargs)
            self.t.append(t[n])  # store corresponding time value
            if n == 0:           # save x once
                np.savez(os.path.join(self.res,'.' + self.filename + '_x.dat'), x=x)        
      
        # Save solution u to a file using numpy.savez
        if self.screen_movie is not False :
    
            # Plot u and mark medium x=x_L and x=x_R
            x_L, x_R = self.medium
            umin, umax = self.yaxis
            title = 'Nx=%d' % (x.size-1)
            if self.title:
                title = self.title + ' ' + title
            if n == 0:
                self.plt.ion()
                self.lines = self.plt.plot(x, u, 'r-',
                    [x_L, x_L], [umin, umax], 'k--',
                    [x_R, x_R], [umin, umax], 'k--')
                self.plt.axis([x[0], x[-1], umin, umax])
                self.plt.xlabel('x')
                self.plt.ylabel('u')
                self.plt.title(title)
                self.plt.legend(['t=%.3f' % t[n]])
                self.plt.pause(0.0001) # Charles
            else:
                # Update new solution
                self.lines[0].set_ydata(u)
                self.plt.legend(['t=%.3f' % t[n]])
                self.plt.draw()
                self.plt.pause(0.0001) # Charles

            # pause
            if t[n] == 0:
                time.sleep(1)  # let initial condition stay 2 s
            else:
                if self.pause is None:
                    pause = 0.0001 if u.size < 100 else 0
                time.sleep(pause)
    

            self.plt.savefig(os.path.join(self.res,'frame_%04d.png' % (n/self.skip)))

    
            if n == (len(t) - 1):   # finished with this run, close plot
                self.plt.close()
        else:
            print(n)