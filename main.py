import numpy as np
import time
N=100

EARTH_RADIUS = 6.371*1e6 #radius in meters
THEIA_RADIUS = 3.05*1e6
EARTH_MASS = 5.972*1e24 #mass in kg
THEIA_MASS = 6.417*1e23

G = 6.673*1e-11 #gravitational constant
REPULSIVE_G = -10*G #used when a sphere is inside another to ensure they don't get stuck
FUSION_DAMPING = 1#Ratio of speed lost when particles fusion
FUSION_POINT= 2000 #If relative speed is inferior to that the 2 sphere fusion
REPULSIVE_ACCEL = -100

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
    speeds = np.zeros([nBalls, 3])
    #positions = []
    #speeds = []
    #for i in range(len(positionsX)):
    #    positions.append(np.array([positionsX[i],positionsY[i], positionsZ[i]]))
    #    speeds.append(np.zeros(3))
    positions[:, 0] = positionsX
    positions[:, 1] = positionsY
    positions[:, 2] = positionsZ
    masses = np.zeros(nBalls)+ballMass

    radiuses = np.zeros(nBalls)+ballRadius

    return positions, masses, speeds, radiuses

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


def computeDeltas(positions, speeds, masses, radiuses, dt):
    P = np.array(positions)
    V = np.array(speeds)
    M = np.array(masses)
    deltaV = np.zeros(V.shape)
    accels = np.zeros(P.shape)
    mask = np.ones(len(P), dtype=np.bool)
    for i in range(len(P)):
        mask[i-1]=1
        mask[i]=0
        P1 = P[i]
        r1 = radiuses[i]
        dirs = (P-P1)[mask]
        redRadiuses = radiuses[mask]
        distances = np.linalg.norm(dirs, axis=1)
        distances = np.where(distances > r1+redRadiuses, distances, r1+redRadiuses)
        gees = G*M[mask]/distances**2
        normals = dirs/distances.reshape(len(dirs), 1)
        accelerations = gees.reshape(len(dirs), 1)*normals
        g = np.sum(accelerations, axis=0)
        deltaV[i] = g*dt
        accels[i] = g
    return deltaV, accels


def computeCollision(positions, speeds, masses, radiuses, i, j):

    pos1 = positions[i]
    pos2 = positions[j]
    speed1 = speeds[i]
    speed2 = speeds[j]
    m1 = masses[i]
    m2 = masses[j]

    check1 = np.dot(speed1, speed2)
    normal1 = (pos2-pos1)/np.linalg.norm(pos2-pos1)
    normal2 = -normal1
    checkSpeed1 = np.dot(speed1, normal1)
    checkSpeed2 = np.dot(speed2, normal2)

    addSpeed1 = checkSpeed1*normal1
    addSpeed2 = checkSpeed2*normal2


    if check1 < 0 and checkSpeed1 < 0:
        return "nothing", 0

    speedDiff = speed2-speed1
    speedNorm = np.linalg.norm(speedDiff)

    if speedNorm < FUSION_POINT:
        r1=radiuses[i]
        r2=radiuses[j]
        newM = m1+m2
        newPos = (m1/(newM))*pos1 + (m2/(newM))*pos2
        newSpeed = (m1/(newM))*speed1 + (m2/(newM))*speed2
        newR = np.sqrt((r1**2 + r2**2))
        positions[i] = newPos; speeds[i] = newSpeed; masses[i] = newM; radiuses[i] = newR
        masses[j] = radiuses[j] = 0
        return "fuse", j

    
    DAMPING_REDUCTION=1e6
    damp = np.max([0.9 ,np.exp(-speedNorm/DAMPING_REDUCTION)])
    
    speed1 -= addSpeed1
    speed1 += addSpeed2
    speed2 += addSpeed1
    speed2 -=addSpeed2
    speed1 = damp*speed1
    speed2 = damp*speed2
    speeds[i] = speed1
    speeds[j] = speed2
    return "collided", 0


def solveCollision(p1, p2, v1, v2, r1, r2, dt):
    #d=p2-p1
    c = deg0Term = (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 +  (p2[2]-p1[2])**2 - (r1+r2)**2
    b = deg1Term = 2*(p2[0]-p1[0])*(v2[0]-v1[0]) + 2*(p2[1]-p1[1])*(v2[1]-v1[1]) + 2*(p2[2]-p1[2])*(v2[2]-v1[2])
    a = deg2Term = (v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + (v2[2]-v1[2])**2

    delta = b**2 - 4*a*c

    if delta<0:
        return None
    sqrtDelta = np.sqrt(delta)
    x1 = (-b -sqrtDelta)/(2*a)
    x2 = (-b +sqrtDelta)/(2*a)

    x1Valid = not (x1 > dt or x1 < 0)
    x2Valid = not (x2 > dt or x2 < 0)
    
    if x1Valid and x2Valid:
        return np.min([x1, x2])
    elif x1Valid and not x2Valid:
        return x1
    elif x2Valid and not x1Valid:
        return x2
    else:
        return None
    

def updateSpeedAndPos(positions, speeds, A, dt):
    speeds += A*dt
    positions += speeds*dt + A * dt**2

def handleUpdate(positions, speeds, newSpeeds, accelerations, masses, radiuses, dt):
    
    collisionsRoots=[]
    collisionsIndexes=[]
    hasBeenTreated = np.identity(len(positions))

    for i in range(len(positions)):
        for j in range(len(positions)):
            if hasBeenTreated[i, j]==1:
                continue
            hasBeenTreated[i, j]=1
            root = solveCollision(positions[i], positions[j], newSpeeds[i], newSpeeds[j], radiuses[i], radiuses[j], dt)
            if root is not None:
                collisionsRoots.append(root)
                collisionsIndexes.append([i, j])

    sortedIndexes = np.argsort(collisionsRoots)
    elapsedTime = 0
    treated = np.zeros(len(positions))
    fused = []

    for i in range(len(collisionsRoots)):

        index = sortedIndexes[i]
        newDt = collisionsRoots[index]-elapsedTime
        firstBallIndex = collisionsIndexes[index][0]
        secondBallIndex = collisionsIndexes[index][1]
        if treated[firstBallIndex] == 1 or treated[secondBallIndex] == 1:
            continue
        treated[firstBallIndex] = 1 
        treated[secondBallIndex] = 1
        updateSpeedAndPos(positions, speeds, accelerations, newDt)
        elapsedTime += newDt
        print(elapsedTime)
        action, index = computeCollision(positions, speeds, masses, radiuses, firstBallIndex, secondBallIndex)
        if action == "fuse":
            fused.append(index)

    updateSpeedAndPos(positions, speeds, accelerations, dt-elapsedTime)


    if len(fused)==0:
        pass
    else:
        mask = np.ones(len(positions), dtype=np.bool)
        for index in fused:
            mask[index]=0
        positions = positions[mask]
        speeds = speeds[mask]
        masses = masses[mask]
        radiuses = radiuses[mask]

    return positions, speeds, masses, radiuses
        

def iteration(positions, speeds, masses, radiuses, dt):
    V = np.array(speeds)
    dV, A = computeDeltas(positions, speeds, masses, radiuses, dt)
    newSpeeds = np.zeros(V.shape)
    newSpeeds = speeds + dV
    return handleUpdate(positions, speeds, newSpeeds, A, masses, radiuses, dt)
                

import tkinter

top = tkinter.Tk()

SIZE_X=800
SIZE_Y=800

C = tkinter.Canvas(top, bg="white", height=SIZE_Y, width=SIZE_X)

UPPER_X = -12*EARTH_RADIUS
UPPER_Y = UPPER_X #Has to be a square
BOTTOM_X = 8*EARTH_RADIUS
BOTTOM_Y = BOTTOM_X

def drawBalls(positions, speeds, radiuses, canvas):

    for i, pos in enumerate(positions):
        
        posX = pos[0]
        posY = pos[1]
        #this is the center of the bawl
        planeX = ( (posX - UPPER_X)/(BOTTOM_X-UPPER_X) )*SIZE_X
        planeY = ( (posY - UPPER_Y)/(BOTTOM_Y-UPPER_Y) )*SIZE_Y
        planeRadius = (radiuses[i] /(BOTTOM_X-UPPER_X))*SIZE_X

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


#drawBalls(positions, speeds, radiuses, C)

C.pack()

#speeds[0, 0]=-5000
#speeds[1, 0]=-3000




earthPositions, earthMasses, earthSpeeds, earthRadiuses = getPlanet(EARTH_RADIUS, 3, EARTH_MASS, 0.74, 0, [0, 0, 0], [0, 0, 0])

nBallsTheia = int(THEIA_MASS/earthMasses[0])
nBallsTheia =1

#earthMasses[0]=earthMasses[1]/2
#earthRadiuses[0]=earthRadiuses[1]/2

theiaPositions, theiaMasses, theiaSpeeds, theiaRadiuses = getPlanet(THEIA_RADIUS, nBallsTheia, THEIA_MASS, 0.74, 0, [-9*EARTH_RADIUS, 0, 0], [0, 0, 0])

counter=0

while counter<100000:
    counter+=1
    drawBalls(earthPositions, earthSpeeds, earthRadiuses, C)
    drawBalls(theiaPositions, theiaSpeeds, theiaRadiuses, C)
    top.update_idletasks()
    top.update()
    earthPositions, earthSpeeds, earthMasses, earthRadiuses = iteration(earthPositions, earthSpeeds, earthMasses,earthRadiuses, 3)
    #iteration(theiaPositions, theiaSpeeds, theiaMasses,theiaRadiuses, 3)
    C.delete("all")

FUSION_DAMPING=0.98


positions=[]
masses=[]
speeds=[]
radiuses=[]

positions.extend(earthPositions); positions.extend(theiaPositions)
masses.extend(earthMasses); masses.extend(theiaMasses)
speeds.extend(earthSpeeds); speeds.extend(theiaSpeeds)
radiuses.extend(earthRadiuses); radiuses.extend(theiaRadiuses)


while 1:
    drawBalls(positions, speeds, radiuses, C)
    top.update_idletasks()
    top.update()
    iteration(positions, speeds, masses, radiuses, 1)
    C.delete("all")



