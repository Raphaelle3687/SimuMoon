import numpy as np
import time
N=100

EARTH_RADIUS = 6.371*1e6 #radius in meters
THEIA_RADIUS = 3.05*1e6
EARTH_MASS = 5.972*1e24 #mass in kg
THEIA_MASS = 6.417*1e23

G = 6.673*1e-11 #gravitational constant
REPULSIVE_G = -10*G #used when a sphere is inside another to ensure they don't get stuck
FUSION_DAMPING = 0.999#Ratio of speed lost when particles fusion
DAMPING_COLLISION_DETECT = 20 #Number of times the acceleration before we use damping
BALL_RADIUS = None
DAMPING_REDUCTION = 1e6
REPULSIVE_ACCEL = -500

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

earthPositions, earthMasses, earthSpeeds, BALL_RADIUS = getPlanet(EARTH_RADIUS, 30, EARTH_MASS, 0.74, 0, [0, 0, 0], [0, 0, 0])

nBallsTheia = int(THEIA_MASS/earthMasses[0])

theiaPositions, theiaMasses, theiaSpeeds, BALL_RADIUS = getPlanet(THEIA_RADIUS, nBallsTheia, THEIA_MASS, 0.74, 0, [-9*EARTH_RADIUS, 0, 0], [0, 0, 0])

positions = earthPositions
masses = earthMasses
speeds = earthSpeeds

positions = appendNPArrayVectors(positions, theiaPositions)
masses = appendNPArray(masses, theiaMasses)
speeds = appendNPArrayVectors(speeds, theiaSpeeds)


def checkCollisions(pos1, speed1, positions, speeds, distances):

    for j, dist in enumerate(distances):
        
        if dist<2*BALL_RADIUS and dist!=0:

    

            #ADD A CHECK TO AVOID DOUBLE REVERSAL/UNNECCESSARY ONES

            pos2 = positions[j]
            speed2 = speeds[j]

            check1 = np.dot(speed1, speed2)
            normal1 = (pos2-pos1)/np.linalg.norm(pos2-pos1)
            normal2 = -normal1
            checkSpeed1 = np.dot(speed1, normal1)
            checkSpeed2 = np.dot(speed2, normal2)

            addSpeed1 = checkSpeed1*normal1
            addSpeed2 = checkSpeed2*normal2


            if check1 < 0 and checkSpeed1 < 0:
                return None, None, None, False

            speedDiff = speed2-speed1
            speedNorm = np.linalg.norm(speedDiff)

            

            damp = np.max([0.9 ,np.exp(-speedNorm/DAMPING_REDUCTION)])
            

            normal = distances[j]/np.linalg.norm(distances[j])
            #correctionVector1 = normal*np.dot(normal, speed1)
            #correctionVector2 = -normal*np.dot(-normal, speed2)
            
            #THIS IS TEMPORARY!! CAN NOT BE PARALELISED IN THIS FORM

            speed1 -= addSpeed1
            speed1 += addSpeed2
            speed2 += addSpeed1
            speed2 -=addSpeed2
            speed1 = damp*speed1
            speed2 = damp*speed2

            return speed1, speed2, j, True
    
    return None, None, None, False

        


def iteration(positions, speeds, masses, dt):

    mask = np.ones(len(positions), dtype=np.bool)
    for i, pos in enumerate(positions):
        mask[i-1]=1
        mask[i]=0
        directions = positions-pos
        distances = np.linalg.norm(directions, axis=1)

        #speed1, speed2, index, hasCollide = checkCollisions(pos, speeds[i], positions, speeds, distances)

        #recoil = np.where(distances>BALL_RADIUS, 1, -REPULSIVE_G)[mask]
        #distances = np.where(distances>BALL_RADIUS, distances, BALL_RADIUS)

        maskedMasses = masses[mask]
        maskedDistances = distances[mask]
        #gees = np.where(maskedDistances > BALL_RADIUS*2, G*((maskedMasses)/maskedDistances**2), REPULSIVE_G*((maskedMasses)/maskedDistances**2) )
        gees = np.where(maskedDistances > BALL_RADIUS*2, G*((maskedMasses)/maskedDistances**2), REPULSIVE_ACCEL )
        damping=1
        if np.any(distances[mask] < BALL_RADIUS*2):
            damping = FUSION_DAMPING

        #gees = G*((masses[mask])/distances[mask]**2)
        directionVector = directions[mask]/distances[mask].reshape(len(directions[mask]), 1)
        directionGees = gees.reshape(len(directions[mask]), 1)*directionVector
        g = np.sum(directionGees, axis=0)

        #(pos1, speed1, positions, speeds, distances, normG)#

        #if hasCollide:
        #    speeds[i] = speed1
        #    speeds[index] = speed2

        speeds[i] += g*dt
        speeds[i]*=damping
        positions[i] += speeds[i]*dt + (dt**2)*g/2

        

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


while 1:
    drawBalls(positions, speeds, C)
    top.update_idletasks()
    top.update()
    time.sleep(0.00001)
    iteration(positions, speeds, masses, 0.3)
    C.delete("all")



