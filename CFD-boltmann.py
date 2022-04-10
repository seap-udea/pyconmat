# Basado en: https://physics.weber.edu/schroeder/fluids/
#########################################
#Requerimientos
#########################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

#########################################
#Propiedades del fluido
#########################################

#Viscodidad cinemática [uh]
nu = 0.01  

#Velocidad inicial del fluido [uh]
U = 0.1

#########################################
#Dimensiones de la región de integración
#########################################

#Ancho y alto de la región de integración [uh]
factor = 3
alto = factor*80
ancho = factor*200

#########################################
#Cantidades derivadas
#########################################

#Constante de relajación [0,2]
w = 1 / (3*nu + 0.5)

#Pesos
w0 = 4.0/9.0
wt = 1.0/9.0
wd  = 1.0/36.0

#########################################
#Inicializa los campos 
#########################################

#Densidades en cada dirección
n0 = w0 * (np.ones((alto,ancho)) - 1.5*U**2)	
nN = wt * (np.ones((alto,ancho)) - 1.5*U**2)
nS = wt * (np.ones((alto,ancho)) - 1.5*U**2)
nE = wt * (np.ones((alto,ancho)) + 3*U + 4.5*U**2 - 1.5*U**2)
nW = wt * (np.ones((alto,ancho)) - 3*U + 4.5*U**2 - 1.5*U**2)
nNE = wd * (np.ones((alto,ancho)) + 3*U + 4.5*U**2 - 1.5*U**2)
nSE = wd * (np.ones((alto,ancho)) + 3*U + 4.5*U**2 - 1.5*U**2)
nNW = wd * (np.ones((alto,ancho)) - 3*U + 4.5*U**2 - 1.5*U**2)
nSW = wd * (np.ones((alto,ancho)) - 3*U + 4.5*U**2 - 1.5*U**2)

#Densidad total
rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW

#Velocidad
ux = (nE + nNE + nSE - nW - nNW - nSW) / rho
uy = (nN + nNE + nNW - nS - nSE - nSW) / rho

#########################################
#Obstáculo
#########################################
#Matriz de la barrera
obstaculo = np.zeros((alto,ancho), bool)

#Barrera lineal
obstaculo[int((alto/2)-factor*8):int((alto/2)+factor*8),int(alto/2)] = True
colorObstaculo=255
if factor==3:
    obstaculo=np.loadtxt("obstaculo.txt").astype("bool")
    colorObstaculo=127

#Crear arreglos booleanos para describir comportamiento de las barreras
obstaculoN = np.roll(obstaculo,  1, axis=0)
obstaculoS = np.roll(obstaculo, -1, axis=0)
obstaculoE = np.roll(obstaculo,  1, axis=1)
obstaculoW = np.roll(obstaculo, -1, axis=1)
obstaculoNE = np.roll(obstaculoN,  1, axis=1)
obstaculoNW = np.roll(obstaculoN, -1, axis=1)
obstaculoSE = np.roll(obstaculoS,  1, axis=1)
obstaculoSW = np.roll(obstaculoS, -1, axis=1)

#########################################
#Colision
#########################################
def colisiones():
    global rho, ux, uy, n0, nN, nS, nE, nW, nNE, nNW, nSE, nSW

    #Densidad
    rho = n0 + nN + nS + nE + nW + nNE + nSE + nNW + nSW

    #Velocidades
    ux = (nE + nNE + nSE - nW - nNW - nSW) / rho
    uy = (nN + nNE + nNW - nS - nSE - nSW) / rho

    #Terminos que se usan con frecuencia
    ux2 = ux * ux 
    uy2 = uy * uy
    u2 = ux2 + uy2
    omu215 = 1 - 1.5*u2
    uxuy = ux * uy

    #Actualiza valores de las densidades
    n0 = (1-w)*n0 + w * w0 * rho * omu215
    nN = (1-w)*nN + w * wt * rho * (omu215 + 3*uy + 4.5*uy2)
    nS = (1-w)*nS + w * wt * rho * (omu215 - 3*uy + 4.5*uy2)
    nE = (1-w)*nE + w * wt * rho * (omu215 + 3*ux + 4.5*ux2)
    nW = (1-w)*nW + w * wt * rho * (omu215 - 3*ux + 4.5*ux2)
    nNE = (1-w)*nNE + w * wd * rho * (omu215 + 3*(ux+uy) + 4.5*(u2+2*uxuy))
    nNW = (1-w)*nNW + w * wd * rho * (omu215 + 3*(-ux+uy) + 4.5*(u2-2*uxuy))
    nSE = (1-w)*nSE + w * wd * rho * (omu215 + 3*(ux-uy) + 4.5*(u2-2*uxuy))
    nSW = (1-w)*nSW + w * wd * rho * (omu215 + 3*(-ux-uy) + 4.5*(u2+2*uxuy))
    
    #Fuerza que el flujo al final sea estacionario 
    nE[:,0] = wt * (1 + 3*U + 4.5*U**2 - 1.5*U**2)
    nW[:,0] = wt * (1 - 3*U + 4.5*U**2 - 1.5*U**2)
    nNE[:,0] = wd * (1 + 3*U + 4.5*U**2 - 1.5*U**2)
    nSE[:,0] = wd * (1 + 3*U + 4.5*U**2 - 1.5*U**2)
    nNW[:,0] = wd * (1 - 3*U + 4.5*U**2 - 1.5*U**2)
    nSW[:,0] = wd * (1 - 3*U + 4.5*U**2 - 1.5*U**2)

#########################################
#Flujo
#########################################
def flujo():
    global nN, nS, nE, nW, nNE, nNW, nSE, nSW

    #Desplazamos los valores de acuerdo a su valor: axis = 0 es norte-sur, axis = 1 es este-oeste
    nN  = np.roll(nN,   1, axis=0)
    nNE = np.roll(nNE,  1, axis=0)
    nNW = np.roll(nNW,  1, axis=0)
    nS  = np.roll(nS,  -1, axis=0)
    nSE = np.roll(nSE, -1, axis=0)
    nSW = np.roll(nSW, -1, axis=0)
    nE  = np.roll(nE,   1, axis=1)
    nNE = np.roll(nNE,  1, axis=1)
    nSE = np.roll(nSE,  1, axis=1)
    nW  = np.roll(nW,  -1, axis=1)
    nNW = np.roll(nNW, -1, axis=1)
    nSW = np.roll(nSW, -1, axis=1)

    #Las partículas rebotan contra los obstáculos
    nN[obstaculoN] = nS[obstaculo]
    nS[obstaculoS] = nN[obstaculo]
    nE[obstaculoE] = nW[obstaculo]
    nW[obstaculoW] = nE[obstaculo]
    nNE[obstaculoNE] = nSW[obstaculo]
    nNW[obstaculoNW] = nSE[obstaculo]
    nSE[obstaculoSE] = nNW[obstaculo]
    nSW[obstaculoSW] = nNE[obstaculo]

#########################################
#Calculo de la vorticidad
#########################################
def vorticidad(ux, uy):
    w=(np.roll(uy,-1,axis=1) - np.roll(uy,1,axis=1)) - (np.roll(ux,-1,axis=0) - np.roll(ux,1,axis=0))
    return w

#########################################
#Animación
#########################################
import matplotlib.pyplot as plt

#Figura que se actualizará con la animación
fig = plt.figure(figsize=(8,3))
fluidoImagen = plt.imshow(vorticidad(ux, uy),
                        origin='lower',
                        norm=plt.Normalize(-.1,.1), 
			cmap="jet",
                        interpolation='none')

#Imagen del obsaculo
obstaculoMatriz = np.zeros((alto, ancho, 4),np.uint8)
obstaculoMatriz[obstaculo,3] = colorObstaculo
obstaculoImagen = plt.imshow(obstaculoMatriz,
                            origin='lower',
                            interpolation='none')

#Animación
def animacion(it):
    for step in range(2):
        flujo()
        colisiones()
    fluidoImagen.set_array(vorticidad(ux, uy))
    return fluidoImagen, obstaculoImagen

anim = animation.FuncAnimation(fig, animacion, interval=1, blit=True)
plt.show()
