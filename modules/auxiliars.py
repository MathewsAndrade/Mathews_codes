#import warnings
import numpy as np
##################################################################################################################################
def estatistica(dado, unidade=None):
    '''
    Calcula os valores minimo, medio, maximo e 
    diferenca entre o maximo e o minimo de um
    dado.
    
    input
    
    dado: numpy array 1D - vetor de dados.
    unidade: string - unidade dos dados.
    
    output
    
    minimo: float - minimo dos dados.
    media: float - media dos dados.
    maximo: float - maximo dos dados.
    variacao: float - diferenca entre o maximo e o
              minimo dos dados.
    '''
    
    assert dado.size > 1, 'O vetor de dados deve ter mais de um elemento'
    
    minimo = np.min(dado)
    maximo = np.max(dado)
    media = np.mean(dado)
    variacao = maximo - minimo
    
    if (unidade != None):
        print ('     min.: %15.5f %s' % (minimo, unidade) )
        print ('    media: %15.5f %s' % (media, unidade) )
        print ('     max.: %15.5f %s' % (maximo, unidade) )
        print ('var. max.: %15.5f %s' % (variacao, unidade) )
    else:
        print ('     min.: %15.5f' % (minimo) )
        print ('    media: %15.5f' % (media) )
        print ('     max.: %15.5f' % (maximo) )
        print ('var. max.: %15.5f' % (variacao) )
        
    return minimo, media, maximo, variacao
##################################################################################################################################
def inv_power_dist(x,y,data,xi,yi,power):
# entradas: 
# x,y,data = 1D arrays com as coordenadas e o dado em cada posicao
# power = real scalar para a potencia a ser utilizada na interpolacao
# xi,yi = 1D arrays com as coordenadas onde se deseja interpolar
# saida:
# dati = dado interpolado   
     
    # verificando qual o menor dos vetores de entrada: 
    n1 = np.size(x)  # dimension of real data
    n2 = np.size(xi) # dimension of interpolated data
    # criacao de arrays:
    peso = np.zeros(n2)
    dati = np.zeros(n2)
    dist = np.zeros( (n2,n1) )
    dummy = 1e-15
    for i in range(n2):
        dist[i,:] = np.sqrt( (x[:] - xi[i] )**2 + ( y[:] - yi[i] )**2 ) # calculo das distancias 
        for j in range(n1):
            if dist[i,j] == 0.:
                dist[i,j] = dummy
    
        peso[i] = np.sum(1.0 / dist[i,:]**power)
        dati[i] = np.sum(data/dist[i,:]**power)
        dati[i] = dati[i]/peso[i]
        
    return dati
##################################################################################################################################
