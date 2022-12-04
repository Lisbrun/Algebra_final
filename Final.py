#Librerias
import numpy as np
from sympy import *
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math
from sympy import *
from PIL import Image
from tkinter import filedialog

#Matrices

#Funciones de entrada y salida.
def leer_matrizE(n,m):
    M = []
    for i in range(n):
        M.append([])
        for j in range(m):
            M[i].append(int(input(f'Fila {i+1}, Columna {j+1}: ')))
    return M
    
def leer_matrizR(n,m):
    M = []
    for i in range(n):
        M.append([])
        for j in range(m):
            M[i].append(float(input(f'Fila {i+1}, Columna {j+1}: ')))
    print("")
    return M

def imprimir_matriz(M):
    n = len(M)
    for i in range(n):
        for j in range(len(M[i])):
            print(M[i][j], end="\t\t")
        print("")

def leer_arreglo(n):
    B = []
    for i in range(n):
        B.append(float(input(f'Ecuación {i+1}: ')))
    return B

def imprimir_arreglo(B):
    for i in B:
        print('\t',i)

def imprimir_arreglo2(B):
    for i in B:
        print(f'\tX{B.index(i)+1} =',i)


#Funciones de Operaciones 

#Suma.
def suma_matriz(A,B):
    C = []
    for i in range(len(A)):
        C.append([])
        for j in range(len(A[i])):
            C[i].append(A[i][j] + B[i][j])
    return C

#Resta.
def resta_matriz(A,B):
    C = []
    for i in range(len(A)):
        C.append([])
        for j in range(len(A[i])):
            C[i].append(A[i][j] - B[i][j])
    return C

#Multiplicacion Escalar.
def productoEscalar(A, E):
    n= len(A)
    for i in range(n):
        for j in range(len(A[i])):
                A[i][j] = A[i][j]*E
    return A

#Multiplicacion Matrices.
def productoMatriz(A, B):
    n= len(A)
    m = len(B[0])
    P = crear_nula(n,m)
    for i in range(len(P)):
        for j in range(len(P[i])):
            for k in range(len(B)):
                P[i][j] += A[i][k]*B[k][j]
    return P

def crear_nula(n,m):
    P = [[0 for j in range(m)] for i in range(n)]
    return  P

#Determinante
def gauss_jordan_det(A):
    n=len(A)
    piv = 0
    aux = 0
    d = 1
    for i in range(n):
        piv = A[i][i]
        d*=piv
        if piv != 0:
            for k in range(n):
                A[i][k] /= piv
            for j in range(n):
                if i != j:
                    aux = A[j][i]
                    for l in range(n):
                        A[j][l] -= aux * A[i][l]
    return d

#Solucion Systema de Ecuaciones.
def gauss_jordan_sol(A,B):
    n = len(A)
    piv = 0
    aux = 0

    for i in range(n):
        A[i].append(B[i])

    for i in range(n):
        piv = A[i][i]
        for k in range(n + 1):
            A[i][k] /= piv
        for j in range(n):
            if i != j:
                aux = A[j][i]
                for k in range(n + 1):
                    A[j][k] -= aux * A[i][k]

    x = []
    for i in range(n):
        x.append(A[i][-1])        
    return x

#Inversa
def identidad(A):
    I = []
    for i in range(len(A)):
        I.append([])
        for j in range(len(A)):
            if i != j:
                I[i].append(0)
            else:
                I[i].append(1)
    return I

def gauss_jordan_inv(A):
    n=len(A)
    piv = 0
    aux = 0
    I = identidad(A)
    for i in range(n):
        piv = A[i][i]
        for k in range(n):
            A[i][k] /= piv
            I[i][k] /= piv
        for j in range(n):
            if i != j:
                aux = A[j][i]
                for k in range(n):
                    A[j][k] -= aux * A[i][k]
                    I[j][k] -= aux * I[i][k]
    return I

#Vectores

#Funciones de entrada y salida.
def a_int():
    return np.array([int(x) for x in input("Digite los valores del vector de enteros separados por un espacio: ").split()])
    
def a_flt():
    return np.array([float(x) for x in input("Digite los valores del vector de reales separados por un espacio: ").split()])

#Norma
def norma(A):
  return np.linalg.norm(A)

#Unitario 
def Vector_U(A):
  return A/norma(A)

def nulo(A):
  b=True
  for i in range(len(A)):
      b= b and A[i]==0
  return b

#Producto punto
def productoP(A,B):
  return A.dot(B)

#Angulo entre vectores
def Angulo(A,B):
   T= np.arccos( np.dot(A,B)/(np.linalg.norm(A)* np.linalg.norm(B) ))
   return np.rad2deg(T)

#Poryeccion
def productoE(A,B):
  for i in range(len(A)):
    A[i]=A[i]*B
  return A

def proyeccion(A,B):
   T= productoE(B,productoP(A,B)/(norma(B))**2)
   return T

#Producto cruz
def Cruz(A,B):
   return np.cross(A,B)


#Main de cada Operacion.

#Matrices
def mainS():
    print("\n")
    print("SUMA".center(35,"·"),"\n")
    n = int(input('Digite el número de filas de las matrices: '))
    m = int(input('Digite el número de columnas de las matrices: '))
    print('\nMatriz A:')
    A = leer_matrizR(n,m)
    imprimir_matriz(A)    
    print('\nMatriz B:')
    B = leer_matrizR(n, m)
    imprimir_matriz(B)
    S = suma_matriz(A, B)
    print('\nMatriz suma:\n')
    imprimir_matriz(S)


def mainR():
    print("\n")
    print("RESTA".center(35,"·"),"\n")
    n = int(input('Digite el número de filas de las matrices: '))
    m = int(input('Digite el número de columnas de las matrices: '))
    print('\nMatriz A:')
    A = leer_matrizR(n,m)
    imprimir_matriz(A)    
    print('\nMatriz B:')
    B = leer_matrizR(n, m)
    imprimir_matriz(B)
    S = resta_matriz(A, B)
    print('\nMatriz suma:\n')
    imprimir_matriz(S)
#mainR()

def mainME():
    print("\n")
    print("MULTIPLICACION CON ESCALAR".center(35, "·"), "\n")
    n = int(input('Digite el número de filas de la matriz: '))
    m = int(input('Digite el número de columnas de la matriz: '))
    print('\nMatriz:')
    A = leer_matrizR(n, m)
    print()
    imprimir_matriz(A)
    E= float(input("\nIngrese el escalar por el cual desea multiplicar la matriz: "))
    print("\n",E,"\n\nMatriz producto:\n")
    imprimir_matriz(productoEscalar(A,E))
#mainME()

def mainMM():
    print("\n")
    print("MULTIPLICACION MATRICES".center(35, "·"), "\n")
    correcto = False
    while correcto == False:
        print('Para multiplicar dos matrices el número de columnas de la matriz A debe ser igual al número de filas de la matriz B')
        n = int(input('Digite el número de filas de la matriz A: '))
        m = int(input('Digite el número de columnas de la matriz A: '))
        p = int(input('Digite el número de filas de la matriz B: '))
        q = int(input('Digite el número de columnas de la matriz B: '))
        if m == p:
            correcto = True
        else:
            print('\nEl número de columnas de A debe ser igual al número de filas de B\n')
    print('\nMatriz A:')
    A = leer_matrizR(n, m)
    print()
    imprimir_matriz(A)
    print('\nMatriz B:')
    B = leer_matrizR(p, q)
    print()
    imprimir_matriz(B)
    print('\nMatriz producto:\n')
    imprimir_matriz(productoMatriz(A,B))
#mainMM()

def mainD():
    print("\n")
    print("DETERMINANTE".center(35, "·"), "\n")
    print('La matriz debe tener igual número de filas y columnas')
    n = int(input('Digite el número de filas de la matriz: '))
    print('\nMatriz:')
    A = leer_matrizR(n,n)
    print()
    imprimir_matriz(A)
    print(f'\nDeterminante de la matriz:',f'{gauss_jordan_det(A):.2f}')
#mainD()

def mainSSE():
    print("\n")
    print("SOLUCION SISTEMA DE ECUACIONES".center(35, "·"), "\n")
    try:
        print('La matriz debe tener igual número de filas y columnas')
        n = int(input('Digite el número de filas de la matriz: '))
        print('\nMatriz:')
        A = leer_matrizR(n, n)
        print()
        imprimir_matriz(A)

        print('\nDigite los resultados del sistema de ecuaciones:\n')
        B = leer_arreglo(n)
        print()
        imprimir_arreglo(B)

        print(f'\nSolución para el sistema de ecuaciones:\n')
        x = gauss_jordan_sol(A,B)
        imprimir_arreglo2(x)

    except ZeroDivisionError:
        print('El sistema no tiene solución')
    except ValueError:
        print('El valor ingresado no es un número')
#mainSSE()

def mainI():
    print("\n")
    print("INVERSA".center(35, "·"), "\n")
    try:
        print('La matriz debe tener igual número de filas y columnas')
        n = int(input('Digite el número de filas de la matriz: '))
        print('\nMatriz:')
        A = leer_matrizR(n, n)
        print()
        imprimir_matriz(A)
        print(f'\nMatriz inversa:\n')
        imprimir_matriz(gauss_jordan_inv(A))
        print("\nComprobacion\n")
        imprimir_matriz(productoMatriz(A,gauss_jordan_inv(A)))
    except ZeroDivisionError:
        print('La matriz es singular. No tiene inversa')
    except ValueError:
        print('El valor ingresado no es un número')
#mainI()

#Vectores

def mainN():
    print("NORMA".center(35,"·"),"\n")
    A = a_flt()
    print(f"La norma de {A} es:", norma(A))
#mainN()  

def mainU():
  print("VECTOR UNITARIO".center(35,"·"),"\n")  
  switch=True
  while switch:
    A = a_flt()
    if nulo(A):
      print("\nAl vector nulo, no se le puede realizar esta operacion")
    else:
      print(f"\nEl vector unitario de  {A} es:",Vector_U(A))    
      print("\nComprobacion:", productoP(Vector_U(A),Vector_U(A)))
      switch=False
#mainU()

def mainPP():
  print("PRODUCTO PUNTO".center(35,"·"),"\n")
  switch= True
  while switch:
    #A = a_int()
    print("Los vectores deben tener la misma cantidad de elmentos.\n\nPrimer vetor")
    A = a_flt()
    print("\nSegundo vector")
    B = a_flt()
    Al=len(A)
    Bl=len(B)
    if (Al==Bl):
      print(f"\nEl producto punto entre {A} y {B} es:",productoP(A,B))
      switch=False
    else:
      print()
#mainPP()

def mainA():
  print("ANGULO ENTRE VECTORSRES".center(35,"·"),"\n")
  switch=True
  while switch:
    print("Los vectores deben tener la misma cantidad de elmentos.\n\nPrimer vetor")
    A = a_flt()
    print("\nSegundo vector")
    B = a_flt()
    Al=len(A)
    Bl=len(B)
    if (Al==Bl):
      print(f"\nEl Angulo punto entre {A} y {B} es: "+"{:.2f}".format(Angulo(A,B))+"°")
      switch=False
    else:
      print()
#mainA()

def mainPro():
  print("PROYECCION".center(35,"·"),"\n")
  switch=True
  while switch:
    print("Los vectores deben tener la misma cantidad de elmentos.\n\nPrimer vetor")
    A = a_flt()
    print("\nSegundo vector")
    B = a_flt()
    Al=len(A)
    Bl=len(B)
    if (Al==Bl):
      print(f"\nLa proyeccion de {A} en {B} es: ", proyeccion(A,B))
      print("\nComprobacion: ",Angulo(proyeccion(A,B),B))
      switch=False
    else:
      print()
#mainPro()

def mainPC():
  print("PRODUCTO CRUZ".center(35,"·"),"\n")
  switch=True
  while switch:
    print("Esta operacion es solo con vectores en R3 y los vectores deben tener la misma cantidad de elmentos.\n\nPrimer vetor")
    A = a_flt()
    print("\nSegundo vector")
    B = a_flt()
    Al=len(A)
    Bl=len(B)
    if (Al==Bl):
      print(f"\nEl producto cruz entre {A} y {B} es:", Cruz(A,B))
      print("\nComprobacion: ", productoP(Cruz(A,B),A))
      switch=False
    else:
      print()


#Valores y vectores propios
def mainVP():
  print("VALORES Y VECTORES PROPIOS".center(35,"·"),"\n")
  print('Ingrese la Matriz la cual desea encontrar valores y vectores propios\n')
  n =int((input('Digite el número de filas de la matriz: ')))
  m = int((input('Digite el número de columnas de la matriz: ')))
  A=leer_matrizR(n,m)
  imp=Matrix(A)
  print("\n")
  imp.eigenvects()

#-----------------------------------Transformaciones y parte visual ----------------------------------------------

def graficar_plano (Polinomio1,Polinomio2,Polinomio3):
    x,y,z = symbols('x y z')
    EcuacionA=Polinomio1[0]*x+Polinomio1[1]*y+Polinomio1[2]*z+Polinomio1[3]
    EcuacionB=Polinomio2[0]*x+Polinomio2[1]*y+Polinomio3[2]*z+Polinomio2[3]
    EcuacionC=Polinomio3[0]*x+Polinomio3[1]*y+Polinomio3[2]*z+Polinomio3[3]


    Solucion = list(linsolve([EcuacionA,EcuacionB,EcuacionC],(x,y,z)))
    x,y = np.linspace(0,10,10),np.linspace(0,10,10)
    vector_cero = [0,0,0]
    X,Y = np.meshgrid(x,y)
    z1 = (Polinomio1[3]+Polinomio1[0]*X+Polinomio1[1]*Y)/Polinomio1[2]
    z2 = (Polinomio2[3]+Polinomio2[0]*X+Polinomio2[1]*Y)/Polinomio2[2]
    z3 = (Polinomio3[3]+Polinomio3[0]*X+Polinomio3[1]*Y)/Polinomio3[2]

    fig = plt.figure()
    ax = fig.add_subplot (111, projection='3d')
    ax.plot_surface(X,Y,z1, alpha= 0.5, cmap = cm.Accent, rstride=10, cstride=10)
    ax.plot_surface(X,Y,z2, alpha= 0.5, cmap = cm.Paired, rstride=10, cstride=10)
    ax.plot_surface(X,Y,z3, alpha= 0.5, cmap = cm.Pastel1, rstride=10, cstride=10)
    ax.set_title("Plano", fontsize=15,fontstyle="oblique",fontweight="bold")
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    plt.show()

def convertir_imagen():
    Imagen = filedialog.askopenfilename(title="Buscar", filetypes=(("Archivos de imagen",".png"), ("Archivos de imagen", ".gif")), initialdir="C:")
    Imagen = Image.open(Imagen)
    Imagen = Imagen.convert("RGB")
    nombre = input("ingrese el nombre con el que desea guardar el archivo:  ")
    Imagen.save(nombre+".jpg")

def Transformacion_linealY(vector,escalar):
    vector = vector
    o = np.array([[0, 0], [0, 0]])

    Matriz_rotacionY=np.array([[1, 0], [0, -1]])
    Resultado= np.dot(Matriz_rotacionY,vector)
    plt.quiver(*o, vector[0],vector[1], color=['green'], scale=15)
    plt.quiver(*o, Resultado[0],Resultado[1], color=['red'], scale=15)


    Matriz_RotacionX= np.array([[-1,0],[0,1]])
    resultado2 = np.dot(Matriz_RotacionX,vector)
    plt.quiver(*o,resultado2[0],resultado2[1], color=['cyan'], scale=15)

    Matriz_escalar= np.array([[float(escalar),0],[0,float(escalar)]])
    resultado3 = np.dot(Matriz_escalar,vector)
    plt.quiver(*o,resultado3[0],resultado3[1], color=['yellowgreen'],scale=15)
    
    plt.quiver(*o,calcular_ortogonal(vector,2)[0],calcular_ortogonal(vector,2)[1], color=['orange'], scale=15)

    leyend= ["Vector Original","Vector Rotacion Y","Vector Rotacion X","Vector Escalar","Vector Ortogonal"]
    plt.legend(loc="upper left",labels=leyend, fontsize=8)
    plt.xlim(-10,10),plt.ylim(-40,40)
    plt.title("Transformaciones Lineales en R2",fontsize=15,fontstyle="oblique",fontweight="bold")
    plt.xlabel("X"), plt.ylabel("Y")
    plt.show()

def calcular_ortogonal(Vector,Y):
    Coord_x= -Vector[1]*Y/Vector[0]
    Coord_y=Y
    print(Coord_x,Coord_y)
    return [Coord_x,Coord_y]


def graficar_3d(vector,escalar,angulo):
    fig = plt.figure()
    ax= fig.add_subplot(projection='3d')
    ax.set_xlim([-1,10])
    ax.set_ylim([-10,10])
    ax.set_zlim([-5,10])
    centro= [0,0,0]
    
    ax.quiver(centro[0],centro[1],centro[2],vector[0],vector[1],vector[2])
    
    matriz_cambio_escala= np.array([[float(escalar),0,0],[0, float(escalar),0],[0,0,1]])
    Resultado= np.dot(matriz_cambio_escala,vector)
    ax.quiver(centro[0],centro[1],centro[2],Resultado[0],Resultado[1],Resultado[2],color="red")


    matriz_Cambio_angulo= np.array([[math.cos(angulo),math.sin(angulo),0],[-math.sin(angulo),math.cos(angulo),0],[0,0,1]])
    Resultado2= np.dot(matriz_Cambio_angulo,vector)
    ax.quiver(centro[0],centro[1],centro[2],Resultado2[0],Resultado2[1],Resultado2[2],color="green",)


    matriz_deformacion =np.array([[1,2,0],[3,1,0],[0,0,1]])
    resultado_deformacion = np.dot(matriz_deformacion,vector)
    ax.quiver(centro[0],centro[1],centro[2],resultado_deformacion[0],resultado_deformacion[1],resultado_deformacion[2],color="yellowgreen", label="Deformacion")

    leyend= ["Vector Original","Vector Escalar","Vector Rotacion","Vector Deformacion"]
    plt.legend(loc="upper left",labels=leyend, fontsize=8)
    ax.set_title("Transformaciones lineales en R3",fontsize=15,fontstyle="oblique",fontweight="bold")
    ax.set_xlabel('X', color='red'), ax.set_ylabel('Y',color='red'), ax.set_zlabel('Z',color='red')
    plt.show()

def seleccionar_imagen():
    imagen = filedialog.askopenfilename(title="Buscar", filetypes=( ("Archivos de imagen", "*.jpg"),("Archivos de imagen", "*.jpeg")), initialdir="C:")
    img= Image.open(imagen)
    matriz_imagen = np.array(img)
    return matriz_imagen



def abrir_imagenes():
    matriz_imagen =seleccionar_imagen()
    grises = escala_grises(matriz_imagen)
    muestra_imagen(grises)
    imagen_volteada=voltear_imagen(matriz_imagen)
    muestra_imagen(imagen_volteada)



def escala_grises(imagen):
    filtro = []
    for i in imagen:
        for r,g,b in i:
            rojo = r*0.2989
            verde = g*0.5870
            azul = b*0.1140
            numero = [int(rojo + verde + azul)]
            filtro.append(numero)

    result = np.reshape(filtro, imagen.shape[:2])
    return result


def muestra_imagen(imagen, cmap="Greys_r"):
    plt.imshow(imagen, cmap=cmap)
    plt.show()


def voltear_imagen(imagen):

    fila = len(imagen[0])
    columna = len(imagen)
    for y in range(columna):
        for x in range (int(fila/2)):
            indice_opuesto = fila - x - 1
            opuesto = imagen[y][indice_opuesto]
            actual = imagen[y][x]
            imagen[y][indice_opuesto] = actual
            imagen[y][x] = opuesto
    muestra_imagen(imagen)

def Transformacion_R2_R3(vector):
    x,y = np.linspace(0,10,10),np.linspace(0,10,10)
    X,Y = np.meshgrid(x,y)
    fig = plt.figure()
    vector_cero=[0,0,0]
    ax = fig.add_subplot (111, projection='3d')

    ax.quiver(vector_cero[0],vector_cero[1],vector_cero[0],vector[0],vector[1],0, color="red")
    Matriz_transformacion= np.array([[1,0],[0,1],[1,1]])
    Transformacion= np.dot(Matriz_transformacion,vector)
    ax.quiver(vector_cero[0],vector_cero[1],vector_cero[0],Transformacion[0],Transformacion[1],Transformacion[2], color="yellowgreen")

    leyend = ["Vector original","Vector transformado"]
    ax.legend(loc="upper left",labels=leyend, fontsize=8)
    ax.set_xlabel('X',color="red"); ax.set_ylabel('Y',color="red"); ax.set_zlabel('Z',color="red")
    ax.set_xlim([-5,5]),ax.set_ylim([-5,5]), ax.set_zlim([-5,5])
    ax.set_title("Transformacion de R2 a R3", fontsize=15,fontstyle="oblique",fontweight="bold")
    plt.show()


def Transformacion_R3_plano(vector):
    x,y = np.linspace(-5,10,2),np.linspace(-5,10,2)
    X,Y= np.meshgrid(x,y)
    fig = plt.figure()   

    z1= (X+Y)-2*vector[0]*X-4*vector[1]*Y/vector[2]
    z2= (X+Y)+2*vector[0]*X-4*vector[1]*Y/vector[2]
    z3= (X+Y)-2*vector[0]*X+4*vector[1]*Y/vector[2]

    ax = fig.add_subplot (111, projection='3d')
    ax.quiver(0,0,0,vector[0],vector[1],vector[2],color='blue', label="vector original")
    ax.plot_surface(X,Y,z1, alpha= 0.5, cmap = cm.Accent, rstride=10, cstride=10)
    ax.plot_surface(X,Y,z2, alpha= 0.5, cmap = cm.Spectral, rstride=10, cstride=10)
    ax.plot_surface(X,Y,z3, alpha= 0.5, rstride=10, cstride=10, cmap= cm.coolwarm)


    leyend = ["Vector original"]
    ax.legend(loc="upper left",labels=leyend, fontsize=8)
    ax.set_xlabel('X',color="red"); ax.set_ylabel('Y',color="red"); ax.set_zlabel('Z',color="red")
    ax.set_title("Transformacion de R3 en plano", fontsize=15,fontstyle="oblique",fontweight="bold")
    plt.show()

def graficar_matriz(matriz):
    fig, ax = plt.subplots()
    ax.matshow(matriz)
    plt.show()


#seleccionar_imagen()
def main():
    iterador= True
    while iterador:
        print("\n------Bienvenido a la parte grafica del curso del Algebra Lineal------\n\n  Seleccione una opcion \n\n 1->Suma de matrices \n 2->Resta de matrices\n 3->Multiplicacion matriz escalar\n 4->Multiplicacion entre matrices \n 5->Determinante  \n 6->Solucion Sistema de ecuaciones-Metodo Guas Jordan \n 7->Matriz inversa\n 8->Norma de vector\n 9->Vector unitario \n 10->Producto punto\n 11->Angulo entre vectores \n 12->Proyeccion\n 13->Producto Cruz \n 14->Graficar en R2 \n 15->Graficar en R3 \n 16->Plano \n 17->Transformacion de R2 a R3 \n 18->Transformacion de R3 a R3 \n 19->Voltear imagen \n 20->Escala de grises \n 21->Graficar matriz \n 22->Valores y vectores propios \n 23->Convertir imagen JPG \n 24->Salir\n")
        seleccion = int(input("seleccion: "))
        
        if seleccion == 1:
            try:
                print()
                mainS()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 2:
            try:
                print()
                mainR()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 3:
            try:
                print()
                mainME()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 4:
            try:
                print()
                mainMM()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 5:
            try:
                print()
                mainD()
            except (ValueError):
                print("--Error en la entrada de datos--")
        
        elif seleccion == 6:
            try:
                print()
                mainSSE()
            except (ValueError):
                print("--Error en la entrada de datos--")
        
        elif seleccion == 7:
            try:
                print()
                mainI()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 8:
            try:
                print()
                mainN()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 9:
            try:
                print()
                mainU()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 10:
            try:
                print()
                mainPP()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 11:
            try:
                print()
                mainA()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 12:
            try:
                print()
                mainPro()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 13:
            try:
                print()
                mainPC()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion == 14:
            try:
                Vector=np.array([int(x) for x in input("Digite los valores del vector en R2 de enteros separados por un espacio: ").split()])
                Escalar= float(input("Digite el escalar: "))
                Transformacion_linealY(Vector,Escalar)
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion==15:
            try:
                Vector=np.array([int(x) for x in input("Digite los valores del vector en R3 de enteros separados por un espacio: ").split()])
                Escalar= float(input("Digite el escalar: "))
                Angulo = float(input("Digite un angulo en radianes: "))
                graficar_3d(Vector,Escalar,Angulo)
            except(ValueError):
                print("--Error en la entrada de datos--")
        
        elif seleccion==16:
            try :
                pol1= np.array([int(x) for x in input("Digite los valores del primer polinomio separados por un espacio: ").split()])
                pol2= np.array([int(x) for x in input("Digite los valores del segundo polinomio separados por un espacio: ").split()])
                pol3= np.array([int(x) for x in input("Digite los valores del tercer polinomio separados por un espacio: ").split()])
                graficar_plano(pol1,pol2,pol3)
            except(ValueError):
                print("--Error en la entrada de datos--")
        elif seleccion==17:
            try:
                Vector=np.array([int(x) for x in input("Digite los valores del vector en R2 de enteros separados por un espacio: ").split()])
                Transformacion_R2_R3(Vector)
            except(ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion==18:
            try:
                Vector=np.array([int(x) for x in input("Digite los valores del vector en R3 de enteros separados por un espacio: ").split()])
                Transformacion_R3_plano(Vector)
            except(ValueError):
                print("--Error en la entrada de datos--")
        elif seleccion==19:
            print("selecciona la imagen")
            matriz =seleccionar_imagen()
            voltear_imagen(matriz)

        elif seleccion==20:
            print("selecciona la imagen")
            matriz =seleccionar_imagen()
            resultado=escala_grises(matriz)
            muestra_imagen(resultado)

        elif seleccion==21:
            matriz=0
            filas = int(input("Digite el numero de filas: "))
            columnas = int(input("Digite el numero de columnas: "))
            matriz = np.empty((filas,columnas),float,order='F')
            print(matriz)
            graficar_matriz(matriz)
            
        elif seleccion == 22:
            try:
                print()
                mainVP()
            except (ValueError):
                print("--Error en la entrada de datos--")

        elif seleccion==23:
            try:
                print()
                convertir_imagen()
            except (ValueError):
                print("--Error en la entrada de datos--")
        elif seleccion==24:
            iterador=False
            print("\nAdios")
        else:
            print("opcion invalida")

if __name__=="__main__":
    main()