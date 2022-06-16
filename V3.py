import numpy as np
from tkinter import *




root = Tk()
root.title('Programm')

root.geometry('630x400')
My_label = Label(root, text = 'Экспертная система',font = ('Times New Roman',16,'bold'))
My_label.grid(row = 0, column = 0)

My_label = Label(root, text = 'Число генераторов:',font = ('Times New Roman',12,'bold'))
My_label.grid(row = 1, column = 0,sticky='W')
numbers = Entry(root, width = 20 )
numbers.grid(row = 1, column = 0 ,sticky='E')
My_label = Label(root, text = 'Номинальные значения частоты генераторов',font = ('Times New Roman',12,'bold'))
My_label.grid(row = 2, column = 0,sticky='W')

numbers2 = Entry(root, width = 20)
numbers2.grid(row = 2, column = 0,sticky='E')
My_label = Label(root, text = 'Номинальная относительная нестабильность частоты генераторов',font = ('Times New Roman',12,'bold'))
My_label.grid(row = 3, column = 0,sticky='W')
numbers3 = Entry(root, width = 20)
numbers3.grid(row = 3, column = 0,sticky='E')
numbers4 = Entry(root, width = 20)
My_label = Label(root, text = 'Количество измерительных интервалов',font = ('Times New Roman',12,'bold'))
My_label.grid(row = 4, column = 0,sticky='W')
numbers4.grid(row = 4, column = 0,sticky='E')
My_label = Label(root, text = 'Номинальная длительность интервала',font = ('Times New Roman',12,'bold'))
My_label.grid(row = 5, column = 0,sticky='W')
numbers5 = Entry(root, width = 20)
numbers5.grid(row = 5, column = 0,sticky='E')






def provvekichina():
    global deltaomageprov
    deltaomageprov = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            deltaomageprov[i][j] = w0 * prov_velichina[i * Kint + j]
    omegaprov = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            omegaprov[i][j] = w0 + deltaomageprov[i][j]
    tprov = []
    for i in range (0,Kint-1):
        tprov.append((t0 * w0)/(w0 + deltaomageprov[i][N-1]))
    deltprov = []
    for i in range (0,Kint-1):
        deltprov = tprov[i]-t0
    global F0prov
    F0prov = w0*t0
    global Fprov
    Fprov = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            Fprov[i][j] = omegaprov[i][j] * tprov[i]



def normraspred():
    global deltaomagenorm
    deltaomagenorm = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            deltaomagenorm[i][j] = w0 * delw0[i * Kint + j]
    omeganorm = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            omeganorm[i][j] = w0 + deltaomagenorm[i][j]
    tnorm = []
    for i in range (0,Kint-1):
        tnorm.append((t0 * w0)/(w0 + deltaomagenorm[i][N-1]))
    deltnorm = []
    for i in range (0,Kint-1):
        deltnorm = tnorm[i]-t0
    global F0
    F0 = w0*t0
    global Fnorm
    Fnorm = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            Fnorm[i][j] = omeganorm[i][j] * tnorm[i]


def ravnraspred():
    global deltaomageravn
    deltaomageravn = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            deltaomageravn[i][j] = w0 * delw01[i * Kint + j]
    omegaravn = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            omegaravn[i][j] = w0 + deltaomageravn[i][j]

    travn = []
    for i in range (0,Kint-1):
        travn.append((t0 * w0)/(w0 + deltaomageravn[i][N-1]))

    deltravn = []
    for i in range (0,Kint-1):
        deltravn = travn[i]-t0

    global Fravn0
    global Fravn
    Fravn0 = w0*t0
    Fravn = [[0  for i in range(N)] for j in range (Kint-1)]
    for i in range (0,Kint-1):
        for j in range (0,N):
            Fravn[i][j] = omegaravn[i][j] * travn[i]

def l0norm():
    delw0norm = [[0  for i in range(R)] for j in range (Kint-1)]
    delw0normlog = []
    for ik in range (0,Kint-1):
        for ir in range (0,R):
            for in1 in range (0,N-1):
                delw0norm[ik][ir] = (Fprov[ik][in1] - F0prov - w0 * deltat[ir])/dt
                L0norm = pow(delw0norm[ik][ir] - deltaomagenorm[ik][in1] ,4)
                
                delw0normlog.append(10 * np.log10(L0norm))
    global minl0norm
    minl0norm = min(delw0normlog)
    print(minl0norm)

def l1norm():
    delw1norm = [[0  for i in range(R)] for j in range (Kint-1)]
    delw1normlog = []
    for ik in range (0,Kint-1):
        for ir in range (0,R):
            for in1 in range (0,N-1):
                delw1norm[ik][ir] = (Fprov[ik][in1] - F0prov - w0 * deltat[ir])/dt
                L1norm = pow(delw1norm[ik][ir] - deltaomageravn[ik][in1] ,4)       
                delw1normlog.append(10 * np.log10(L1norm))
    global minl1norm
    minl1norm = min(delw1normlog)  
    print(minl1norm)
    


def clicked():
    global N,w0,sigma0,Kint,t0,R,dt,delt,delw0,delw01,prov_velichina,deltat

    N  = int(numbers.get())
    w0 = int(numbers2.get())
    sigma0 = float(numbers3.get())
    Kint = int(numbers4.get())
    t0 = float(numbers5.get())
    #N = 1000 #Число генераторов
    #w0 = 2*np.pi*pow(10,9) #номинальные значения частоты генераторов
    #sigma0 = pow(10, -7) #номинальная относительная нестабильность частоты генераторов
    #Kint = 100 #количество измерительных интервалов
    #t0 = 10^(-3) #номинальная длительность интервала


    R = 1000
    dt = (2*np.sqrt(3)*sigma0)/R
    delt = []
    for i in range (0,R):
        delt.append(np.sqrt(3)*sigma0 + dt * i)
    delw0 =  np.random.normal(loc = 0.0, scale = sigma0,size = (N+1)*Kint)
    delw01 = np.random.uniform(-np.sqrt(3)*sigma0,  np.sqrt(3)*sigma0,(N+1)*Kint)
    prov_velichina = np.random.uniform(-np.sqrt(3)*sigma0,  np.sqrt(3)*sigma0,(N+1)*Kint)
    provvekichina()
    normraspred()
    ravnraspred()
    
    dt = (2 * np.sqrt(3) * sigma0)/(10*R)
    deltat = []
    for i in range (0,R):
        deltat.append(((-np.sqrt(3) * sigma0)/10)+dt * i)
    l0norm()
    l1norm()

btn = Button(root, text="Запустить работу программы", command=clicked) 
btn.grid( sticky='EW')
text = Text(width = 50, height = 10,padx = 110) 

def analiz():
    if (minl0norm < minl1norm):
        text.insert("insert","Нормальный закон распределения")
    else:
        text.insert("insert","Равномерный закон распределения")



btn1 = Button(root, text="Анализ", command=analiz) 
btn1.grid( sticky='EW')
text.grid() 
root.mainloop()

