#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from math import log

y_f1 = []
y_f1_promedio = []

y_hitrate = []
y_hitrate_promedio = []

y_recalls = []
y_recalls_promedio_iteracion = []
y_recalls_promedio_clase = []

y_precision = []
y_precision_promedio_iteracion = []
y_precision_promedio_clase = []

y_tiempo = []

# MODIFICAR EN CASO DE QUERER DISTINTOS VALORES EN EL EJE X DEL GRAFICOS
# x_clase sería las iteraciones que se hacen variando el valor que se varíe

x_clase = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100] # Por ejemplo aca varío los vecinos de 10 a 100

K = 10 # K del k-cross-folding
cant_iteraciones = len(x_clase) # modificar x. RANGO DE VARIABLE

# x_iteracion es para gráficar, por cada iteración de la variable siendo modificada, cada iteración del K-cross-folding
# Por eso se calcula como cant_iteraciones * K

x_iteracion = []
a = 10

for j in range(0, cant_iteraciones):
  for i in range (0, K):
    x_iteracion.append(a)

  a += 10

# MODIFICAR EL NOMBRE DE LA CARPETA DE DONDE LEER LOS DATOS (y la extension del archivo, yo flashee y les puse .tex)
# Ejemplo:
# f = open('../resultados/MI_CARPETA/f1.tex', 'r')
promedio = 0
f = open('../resultados/PLSDA_vecinos_10_a_100_dimension_10/f1.tex', 'r')
for j in range(0, cant_iteraciones):
  for i in range(0, K):
    aux = float(f.readline()[:-1])
    promedio += aux
    y_f1.append(aux)

  promedio /= K
  y_f1_promedio.append(promedio)
  promedio = 0
          
f = open('../resultados/PLSDA_vecinos_10_a_100_dimension_10/hitrate.tex', 'r')
promedio = 0
for j in range(0, cant_iteraciones):
  for i in range(0, K):
    aux = float(f.readline()[:-1])
    promedio += aux
    y_hitrate.append(aux)

  promedio /= K
  y_hitrate_promedio.append(promedio)
  promedio = 0

f = open('../resultados/PLSDA_vecinos_10_a_100_dimension_10/recalls.tex', 'r')
prom = 0
for j in range(0, cant_iteraciones):
  for i in range(0, K):
    for z in range(0, 10):
      y_recalls.append(float(f.readline()[:-1]))

    promedio = float(f.readline()[:-1])
    prom += promedio
    y_recalls_promedio_iteracion.append(promedio)

  prom /= K
  y_recalls_promedio_clase.append(prom)
  prom = 0

f = open('../resultados/PLSDA_vecinos_10_a_100_dimension_10/precisiones.tex', 'r')
for j in range(0, cant_iteraciones):
  for i in range(0, K):
    for z in range(0, 10):
      y_precision.append(float(f.readline()[:-1]))

    promedio = float(f.readline()[:-1])
    prom += promedio
    y_precision_promedio_iteracion.append(promedio)

  prom /= K
  y_precision_promedio_clase.append(prom)
  prom = 0


# DESCOMENTAR EN CASO DE HABER HECHO ADEMAS TOMA DE TIEMPOS

#promedio = 0
#f = open('../resultados/PLSDA_vecinos_10_a_100_dimension_10/tiempos.tex', 'r')
#for j in range(0, cant_iteraciones):
#  for i in range(0, K):
#    promedio += float(f.readline()[:-1])
#
#  promedio /= K
#  y_tiempo.append(promedio)
#  promedio = 0

# PARA QUE SE GUARDEN LAS IMAGENES ES NECESARIO CREAR UNA CARPETA DENTRO DE SUS RESULTADOS QUE SE LLAME "graficos"
# Ejemplo:
# plt.savefig('../resultados/MI_CARPETA/graficos/precisiones_clase.png')

# PRIMER IMAGEN

error_config = {'ecolor': '0.3'}
bar_width = 0.96
opacity = 0.4
plt.plot(x_clase, y_precision_promedio_clase, 'ro', #bar_width,
          alpha=opacity,
          color='b',
          #error_kw=error_config,
          label=u"Precisión del categorizador")

plt.xlabel(u"cantidad de vecinos")
plt.ylabel(u"Precisión")
plt.legend()
plt.savefig('../resultados/PLSDA_vecinos_10_a_100_dimension_10/graficos/precisiones_clase.png')
plt.show()

# SEGUNDA IMAGEN

plt.plot(x_iteracion, y_precision_promedio_iteracion, 'ro', #bar_width,
          alpha=opacity,
          color='b',
          label=u"Precisión del categorizador por clase")

plt.xlabel(u"cantidad de vecinos")
plt.ylabel(u"Precisión")
plt.legend()

plt.tight_layout()
plt.savefig('../resultados/PLSDA_vecinos_10_a_100_dimension_10/graficos/precisiones_iteracion.png')
plt.show()

# TERCER IMAGEN

plt.plot(x_clase, y_recalls_promedio_clase, 'ro', #bar_width,
          alpha=opacity,
          color='b',
          label=u"recall del categorizador")

plt.xlabel(u"cantidad de vecinos")
plt.ylabel(u"Recall")
plt.legend()

plt.tight_layout()
plt.savefig('../resultados/PLSDA_vecinos_10_a_100_dimension_10/graficos/recall_clase.png')
plt.show()

# CUARTA IMAGEN

plt.plot(x_iteracion, y_recalls_promedio_iteracion, 'ro', #bar_width,
          alpha=opacity,
          color='b',
          label=u"recall del categorizador por clase")

plt.xlabel(u"cantidad de vecinos")
plt.ylabel(u"Recall")
plt.legend()

plt.tight_layout()
plt.savefig('../resultados/PLSDA_vecinos_10_a_100_dimension_10/graficos/recall_iteracion.png')
plt.show()

# QUINTA IMAGEN

plt.plot(x_clase, y_hitrate_promedio, 'ro', #bar_width,
          alpha=opacity,
          color='b',
          label=u"Hit rate del categorizador")

plt.xlabel(u"cantidad de vecinos")
plt.ylabel(u"Hit rate")
plt.legend()

plt.tight_layout()
plt.savefig('../resultados/PLSDA_vecinos_10_a_100_dimension_10/graficos/hitrate_promedio.png')
plt.show()

# SEXTA IMAGEN

plt.plot(x_clase, y_f1_promedio, 'ro', #bar_width,
          alpha=opacity,
          color='b',
          label=u"F1 del categorizador")

plt.xlabel(u"cantidad de vecinos")
plt.ylabel(u"F1")
plt.legend()

plt.tight_layout()
plt.savefig('../resultados/PLSDA_vecinos_10_a_100_dimension_10/graficos/F1_promedio.png')
plt.show()

# SEPTIMA IMAGEN
# DESCOMENTAR EN CASO DE HACER TOMA DE TIEMPOS

#plt.plot(x_clase, y_tiempo, 'ro', #bar_width,
#          alpha=opacity,
#          color='b',
#          label=u"Tiempo del PLSDA+kNN")
#
#plt.xlabel(u"cantidad de vecinos")
#plt.ylabel(u"tiempo")
#plt.legend()
#
#plt.tight_layout()
#plt.savefig('../resultados/PLSDA_vecinos_10_a_100_dimension_10/graficos/tiempos.png')
#plt.show()

