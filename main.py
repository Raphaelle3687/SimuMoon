import numpy as np
import time
N=100

EARTH_RADIUS = 6.371*1e6 #radius in meters
THEIA_RADIUS = 3.05*1e6
EARTH_MASS = 5.972*1e24 #mass in kg
THEIA_MASS = 6.417*1e24

G = 6.673*1e-11 #gravitational constant
REPULSIVE_G = -10*G #used when a sphere is inside another to ensure they don't get stuck
FUSION_DAMPING = 10#Ratio of speed lost when particles fusion
DAMPING_COLLISION_DETECT = 3500 #Number of times the acceleration before we use damping
BALL_RADIUS = None
DAMPING_REDUCTION = 1e6
REPULSIVE_ACCEL = -100

positions=[]
masses=[]
speeds=[]

def getPlanet(planetRadius, nBalls, planetMass, packingConstant, rotationSpeed, centerPoint, directionSpeed, ballRadius=None):
    
    if ballRadius is None:
        planetVolume = (np.pi)*planetRadius**2
        ballVolume = planetVolume*packingConstant/nBalls
        ballRadius = np.sqrt(ballVolume/np.pi)
    ballMass = planetMass / nBalls

    limits = 2.5*planetRadius
    positionsX = np.random.uniform(low=-limits+centerPoint[0], high=limits+centerPoint[0], size=[nBalls])
    positionsY = np.random.uniform(low=-limits+centerPoint[1], high=limits+centerPoint[1], size=[nBalls])
    positionsZ = np.random.uniform(low=-limits+centerPoint[2], high=limits+centerPoint[2], size=[nBalls])
    positions = np.zeros([nBalls, 3])
    positions[:, 0] = positionsX
    positions[:, 1] = positionsY
    positions[:, 2] = positionsZ
    masses = np.zeros(nBalls)+ballMass
    speeds = np.zeros([nBalls, 3])

    return positions, masses, speeds, ballRadius

def appendNPArray(arr1, arr2):
    newArray = np.zeros(len(arr1)+len(arr2))
    newArray[:len(arr1)]=arr1
    newArray[len(arr1):]=arr2
    return newArray

def appendNPArrayVectors(arr1, arr2):
    newArray = np.zeros([len(arr1)+len(arr2), 3])
    newArray[:len(arr1), :]=arr1
    newArray[len(arr1):, :]=arr2
    return newArray

def pairInteraction(pos1, pos2, speed1, speed2, m1, m2, dt):

    distance = pos2-pos1
    normDistance = np.linalg.norm(distance)
    normal1 = (distance)/normDistance
    damping = 1
    relativeSpeedNorm = np.linalg.norm(speed2-speed1)
    if normDistance > 2*BALL_RADIUS:
        g2 = G*m1/(normDistance**2)
        g1=g2*(m2/m1)
    else:
        g1=REPULSIVE_ACCEL
        g2=g1
        if relativeSpeedNorm > DAMPING_COLLISION_DETECT:
            damping=FUSION_DAMPING
    
    speed1 += dt*g1*(normal1)
    speed2 += dt*g2*(-normal1)
    print(damping)
    speed1*=damping
    speed2*=damping
    pos1 += speed1*dt
    pos2 += speed2*dt
    return pos1, pos2, speed1, speed2


def iteration(positions, speeds, masses, dt):

    mask=np.ones(len(positions))

    for i in range(len(positions)):
        mask[i]=0
        for j in range(len(positions)):
            if mask[j]!=0:

                positions[i], positions[j], speeds[i], speeds[j] = \
                    pairInteraction(positions[i], positions[j], speeds[i], speeds[j], masses[i], masses[j], dt)
                

        

#for i in range(0, 2000):
#    distance = np.sqrt(np.sum((positions[0]-positions[1])**2))
#    print(distance, BALL_RADIUS*2, BALL_RADIUS)
#    iteration(positions, speeds, masses, 5)

import tkinter

top = tkinter.Tk()

SIZE_X=800
SIZE_Y=800

C = tkinter.Canvas(top, bg="white", height=SIZE_Y, width=SIZE_X)

UPPER_X = -12*EARTH_RADIUS
UPPER_Y = UPPER_X #Has to be a square
BOTTOM_X = 8*EARTH_RADIUS
BOTTOM_Y = BOTTOM_X

def drawBalls(positions, speeds, canvas):

    for i, pos in enumerate(positions):
        
        posX = pos[0]
        posY = pos[1]
        #this is the center of the bawl
        planeX = ( (posX - UPPER_X)/(BOTTOM_X-UPPER_X) )*SIZE_X
        planeY = ( (posY - UPPER_Y)/(BOTTOM_Y-UPPER_Y) )*SIZE_Y
        planeRadius = (BALL_RADIUS /(BOTTOM_X-UPPER_X))*SIZE_X

        x1 = planeX-planeRadius
        x2 = planeX+planeRadius
        y1 = planeY-planeRadius
        y2 = planeY+planeRadius

        if i%3==0:
            color="red"
        elif i%3==1:
            color="blue"
        else:
            color="green"

        xLine=speeds[i][0]/(EARTH_RADIUS*1e-5)
        yLine=speeds[i][1]/(EARTH_RADIUS*1e-5)


        canvas.create_oval(x1, y1, x2, y2, fill='grey')

        #canvas.create_line(planeX, planeY, planeX+xLine, planeY+yLine, fill=color, width=5)


drawBalls(positions, speeds, C)

C.pack()

#speeds[0, 0]=-5000
#speeds[1, 0]=-3000


earthPositions, earthMasses, earthSpeeds, BALL_RADIUS = getPlanet(EARTH_RADIUS, 10, EARTH_MASS, 0.74, 0, [0, 0, 0], [0, 0, 0])

nBallsTheia = int(THEIA_MASS/earthMasses[0])

theiaPositions, theiaMasses, theiaSpeeds, BALL_RADIUS = getPlanet(THEIA_RADIUS, nBallsTheia, THEIA_MASS, 0.74, 0, [-9*EARTH_RADIUS, 0, 0], [0, 0, 0])

FUSION_DAMPING=0.02

counter=0

while counter<1000:
    counter+=1
    drawBalls(earthPositions, earthSpeeds, C)
    drawBalls(theiaPositions, theiaSpeeds, C)
    top.update_idletasks()
    top.update()
    iteration(earthPositions, earthSpeeds, earthMasses, 1)
    iteration(theiaPositions, theiaSpeeds, theiaMasses, 3)
    C.delete("all")

FUSION_DAMPING=0.98

positions = earthPositions
masses = earthMasses
speeds = earthSpeeds

positions = appendNPArrayVectors(positions, theiaPositions)
masses = appendNPArray(masses, theiaMasses)
speeds = appendNPArrayVectors(speeds, theiaSpeeds)

while 1:
    drawBalls(positions, speeds, C)
    top.update_idletasks()
    top.update()
    iteration(positions, speeds, masses, 0.1)
    C.delete("all")



