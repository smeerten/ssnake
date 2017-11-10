#!/usr/bin/env python

# Copyright 2016 - 2017 Bas van Meerten and Wouter Franssen

# This file is part of ssNake.
#
# ssNake is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ssNake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ssNake. If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import random
import numpy as np
try:
    from PyQt4 import QtGui, QtCore
    from PyQt4 import QtGui as QtWidgets
    QT = 4
except ImportError:
    from PyQt5 import QtGui, QtCore, QtWidgets
    QT = 5
import matplotlib
if QT ==4:
    matplotlib.use('Qt4Agg')
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
else:
    matplotlib.use('Qt5Agg')
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

class runSnake: #The actual game
    def __init__(self, root):
        self.father = root
        self.dims = self.father.dims
        self.running = False
        self.score = 0
        #self.start()
        self.father.PlotFrame.message()
        self.timer = QtCore.QTimer()
        self.timer.setInterval(150)
        self.timer.timeout.connect(self.moveSnake)

    def start(self):
        self.score = 0
        if self.running: #If a restart, reset the timer
            self.timer.stop()
        self.snake = [np.array([5,5]),np.array([5,4]),np.array([5,3])] #List with starting points of the snake
        self.dir = np.array([0,-1]) #row and column heading. Start with going left
        self.food = None 
        self.placeFood() #Generate a random food position
        self.plot() #Plot the first frame
        self.mainLoop() #Start the game loop
        self.running = True

    def clear(self):
        if self.running: #If a restart, reset the timer
            self.timer.stop()
            self.running  = False
        self.snake = []
        self.food = []
        self.plot()
        self.father.PlotFrame.message(self.score)

    def mainLoop(self):
        self.timer.start()

    def placeFood(self):
        rows = self.dims[0]
        cols = self.dims[1]
        while True:
            row = random.randint(0,rows-1)
            col = random.randint(0,cols-1)
            New = True
            for elem in self.snake:
                if elem[0] == row and elem[1] == col:
                    New = False
            if New:
                break
        self.food = [np.array([row,col])]

    def moveSnake(self):
        NewPos = self.snake[-1] + self.dir 
        #Check hitting wall
        if NewPos[0] == self.dims[0] or NewPos[1] == self.dims[1] or NewPos[0] < 0 or NewPos[1] <0:
            self.clear()
            return
        #Check hitting itself. Ignore the last point of the tail (it just slips away)
        for pos in self.snake[1::]:
            if NewPos[0] == pos[0] and NewPos[1] == pos[1]:
                self.clear()
                return

        if np.allclose(self.food,NewPos):#eatfood
            self.score += 1
            self.snake.append(NewPos)
            self.placeFood()
            self.plot('eat')
        else: #Normal move
            self.snake.append(NewPos)
            self.snake.pop(0)
            self.plot('move')

    def keyPressed(self,event):
        if event.key == "q":
            self.clear()
            self.father.PlotFrame.message()
        if event.key == "r":
            self.start()
            return
        # now process keys that only work if the game is not over
        if event.key == "up":
            self.dir = np.array([+1, 0])
        elif event.key == "down":
            self.dir = np.array([-1, 0])
        elif event.key == "left":
            self.dir = np.array([0, -1])
        elif event.key == "right":
            self.dir = np.array([0, +1])

    def plot(self,type = 'init'):
        self.father.PlotFrame.update(self.snake,self.food,type)
 
class PlotFrame(object):

    def __init__(self, root, fig, canvas):
        self.root = root
        self.fig = fig
        self.canvas = canvas
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)
        self.snakeColours = [np.array([0.0,0.9,0.7,1.0]),np.array([0.0,0.0,1.0,1.0])]

        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()
        #self.ax.axis('equal')
        self.ax.set_xticks([]) 
        self.ax.set_yticks([]) 
        self.ax.set_xlim(-1, self.root.dims[1])
        self.ax.set_ylim(-1, self.root.dims[0] )
        self.canvas.draw()
        self.snakeArtists = []
        self.foodArtists = []
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.text = []

    def message(self,score = 0):
        for elem in self.text:
            elem.remove()
        self.text = []
        if score != 0:
            t2 = 'Game over. Score: ' + str(score) + '\n\n'
        else:
            t2 = ''
        t = "Welcome to Snake in ssNake!\nPress <r> to start and <q> ro return here."
        self.text.append(self.ax.text(15,10,t2 + t ,ha='center',va='center'))
        self.canvas.draw()
  
    def update(self,snake,food,type = 'init'):
        #Clear artist list if nesesarry
        if type == 'init': #Clear artist all
            for elem in self.snakeArtists:
                elem.remove()
            self.snakeArtists = []
            for elem in self.foodArtists:
                elem.remove()
            self.foodArtists = []
        elif type == 'move': #pop first
            self.snakeArtists[0].remove()
            self.snakeArtists.pop(0)
        elif type == 'eat':
            self.foodArtists[0].remove()
            self.foodArtists.pop(0)
            
        #Create new artists (i.e. circles)
        if type == 'init':
            for elem in snake:
                self.snakeArtists.append(self.ax.add_artist(plt.Circle((elem[1], elem[0]), 0.5, color='r')))
            if food is not None:
                for pos in food:
                    self.foodArtists.append(self.ax.add_artist(plt.Circle((pos[1], pos[0]), 0.5, color='g')))
        elif type == 'move':
                self.snakeArtists.append(self.ax.add_artist(plt.Circle((snake[-1][1], snake[-1][0]), 0.5, color='r')))
        elif type == 'eat':
                self.snakeArtists.append(self.ax.add_artist(plt.Circle((snake[-1][1], snake[-1][0]), 0.5, color='r')))
                self.foodArtists.append(self.ax.add_artist(plt.Circle((food[-1][1], food[-1][0]), 0.5, color='g')))
        
        #Set the colors for the snake
        if len(self.snakeArtists) > 0:
            Totlen = len(self.snakeArtists)
            step = (self.snakeColours[1] - self.snakeColours[0]) / (Totlen )
            for art in range(len(self.snakeArtists)):
                self.snakeArtists[art].set_color(tuple(self.snakeColours[0] + step * art))
        #Restore default plot, and add the artists
        self.canvas.restore_region(self.background) 
        for art in self.snakeArtists:
            self.ax.draw_artist(art)
        for art in self.foodArtists:
            self.ax.draw_artist(art)
        self.canvas.blit(self.ax.bbox)

class MainProgram(QtWidgets.QMainWindow):

    def __init__(self, root):
        super(MainProgram, self).__init__()
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.main_widget = QtWidgets.QWidget(self)
        self.mainFrame = QtWidgets.QGridLayout(self.main_widget)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.gca()
        self.mainFrame.addWidget(self.canvas, 0, 0)
        self.mainFrame.setColumnStretch(0, 1)
        self.mainFrame.setRowStretch(0, 1)
        self.dims = [20,30] #Rows,columns
        self.PlotFrame = PlotFrame(self, self.fig, self.canvas)
        self.game = runSnake(self)
        self.canvas.mpl_connect('key_press_event', self.game.keyPressed)

    def fileQuit(self):
        self.close()

if __name__ == '__main__':
    root = QtWidgets.QApplication(sys.argv)
    mainProgram = MainProgram(root)
    root.setWindowIcon(QtGui.QIcon(os.path.dirname(os.path.realpath(__file__)) + '/logo.gif'))
    mainProgram.setWindowTitle(u"Snake in ssNake!")
    mainProgram.show()
    sys.exit(root.exec_())


