#!/home/samhaug/anaconda2/bin/python

from matplotlib import pylab as p
import numpy as np

class DragHandler(object):
    """ A simple class to handle Drag n Drop.

    This is a simple example, which works for Text objects only.
    """
    def __init__(self, figure=None) :
        """ Create a new drag handler and connect it to the figure's event system.
        If the figure handler is not given, the current figure is used instead
        """
        if figure is None : figure = p.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None

        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)

    def on_pick_event(self, event):
        " Store which text object was picked and were the pick event occurs."

        if event.artist.get_label() == 'sim':
            self.dragged = event.artist
            self.xdata = event.artist.get_data()[0]
            self.ydata = event.artist.get_data()[1]
            self.pick_pos = event.mouseevent.xdata
        return True

    def on_release_event(self, event):
        " Update text position and redraw"

        if self.dragged is not None :
            newx = event.xdata
            newy = np.roll(self.ydata,int(newx-self.pick_pos))
            self.dragged.set_data(self.xdata,newy)
            self.dragged = None
            p.draw()
        return True

