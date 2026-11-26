# Visualización 3D de una arteria coronaria en CTA

Herramienta en C++ con VTK para explorar una arteria coronaria a partir de una angiografía por tomografía computarizada (CTA). Extrae la superficie del vaso, la aísla del resto del tejido, estima un radio a lo largo de la línea central y colorea la zona de mayor estrechamiento estimado.

Proyecto final de Computación Gráfica, Pontificia Universidad Javeriana.

> **Esto no es una herramienta de diagnóstico.** El estrechamiento es un indicador visual exploratorio, sin validación contra una referencia clínica.

![Arteria con contexto del corazón](https://github.com/user-attachments/assets/d2d3ffc1-39d4-49f9-90e6-38abd201b808)

## Qué hace

1. Lee un volumen CTA (`.mhd`) y la línea central de un vaso (`.vtk`).
2. Verifica que ambos comparten sistema de coordenadas.
3. Calcula un umbral automático a partir de la intensidad de la arteria.
4. Extrae la superficie con Marching Cubes sobre una caja recortada.
5. Aísla la arteria por distancia a la línea central y por conectividad.
6. Estima el radio en cada punto de la línea y calcula el estrechamiento.
7. Dibuja un tubo sobre la línea, coloreado con el mapa viridis.
8. Opcionalmente dibuja el contexto del corazón y las calcificaciones.

## Requisitos

- C++11 o superior
- CMake 3.12 o superior
- VTK 9.7
- Probado en Ubuntu 24.04 

## Compilación

```bash
mkdir build && cd build
cmake ..
make
```

Ejemplo mínimo de `CMakeLists.txt`.

```cmake
cmake_minimum_required(VERSION 3.12)
project(coronarias)

find_package(VTK REQUIRED)

add_executable(coronarias main.cpp)
target_link_libraries(coronarias ${VTK_LIBRARIES})
vtk_module_autoinit(TARGETS coronarias MODULES ${VTK_LIBRARIES})
```

## Datos

Los datos no están en este repositorio. Provienen del marco de evaluación de Kirişli et al. (2013), que reúne 48 estudios CTA de tres centros y tres fabricantes. 

Cada caso incluye el volumen (`image.mhd` e `image.raw`) y líneas centrales por vaso. Este proyecto usa la del grupo LKEB corregida manualmente (`Team_LKEB_ManuallyCorrected_vessel00.vtk`) en `dataset00` y `dataset08`.

## Uso

```bash
./coronarias image.mhd linea_central.vtk [umbral] [dist_mm] [umbral_corazon] [umbral_calcio]
```

| Argumento | Descripción |
|---|---|
| `image.mhd` | Volumen CTA. Obligatorio. |
| `linea_central.vtk` | Línea central del vaso. Obligatorio. |
| `umbral` | Valor de la isosuperficie. Si se omite o es 0, se usa la mitad de la intensidad en el punto medio de la línea. |
| `dist_mm` | Distancia de recorte a la línea central. Si es mayor que 0, aísla la arteria, calcula el radio y colorea el tubo. Si se omite o es 0, muestra el control inicial. |
| `umbral_corazon` | Si es mayor que 0 (y `dist_mm` también), dibuja el contexto del corazón con ese umbral. |
| `umbral_calcio` | Si es mayor que 0 (y `dist_mm` también), dibuja las calcificaciones cerca de la línea. Debe ser claramente mayor que la intensidad de la arteria. |

### Modos de ejecución

```bash
D=/ruta/a/1_training/dataset08

# Control inicial, línea amarilla sobre la superficie cruda
./coronarias $D/image.mhd $D/Team_LKEB_ManuallyCorrected_vessel00.vtk

# Arteria aislada con tubo coloreado
./coronarias $D/image.mhd $D/Team_LKEB_ManuallyCorrected_vessel00.vtk 0 6

# Arteria dentro del corazón (lento, del orden de 2 millones de triángulos)
./coronarias $D/image.mhd $D/Team_LKEB_ManuallyCorrected_vessel00.vtk 0 6 250

# Arteria con calcificaciones (el umbral hay que calibrarlo por caso)
./coronarias $D/image.mhd $D/Team_LKEB_ManuallyCorrected_vessel00.vtk 0 6 0 600
```

### Controles

Estilo `vtkInteractorStyleTrackballCamera`.

| Acción | Control |
|---|---|
| Rotar | Botón izquierdo y arrastrar |
| Desplazar | Botón central, o Mayús y botón izquierdo |
| Zoom | Botón derecho y arrastrar, o rueda |
| Salir | Tecla `q` |

## Cómo funciona

```
Lectura -> Umbral -> Recorte de caja -> Marching Cubes
        -> Distancia a la línea -> Clip -> Conectividad
        -> Radio -> Suavizado -> Percentil 90 -> Estrechamiento
        -> Tubo coloreado -> Escena
```

| Etapa | Clases de VTK |
|---|---|
| Lectura | `vtkMetaImageReader`, `vtkPolyDataReader` |
| Recorte del volumen | `vtkExtractVOI` |
| Superficie | `vtkMarchingCubes` |
| Distancia a la línea | `vtkPointLocator` |
| Aislamiento | `vtkClipPolyData`, `vtkPolyDataConnectivityFilter` |
| Radio | `vtkPointLocator` sobre la pared aislada |
| Color | `vtkTubeFilter`, `vtkColorTransferFunction` |
| Escena | `vtkPolyDataMapper`, `vtkActor`, `vtkRenderer`, `vtkRenderWindow` |

### Umbral

La intensidad de la sangre con contraste cambia entre casos (324 en `dataset00`, 345 en `dataset08`), así que un valor fijo no sirve. El programa lee la intensidad en el punto medio de la línea y usa la mitad.

### Estrechamiento

El radio en cada punto de la línea es la distancia al vértice más cercano de la superficie aislada. Se suaviza promediando cada valor con sus dos vecinos a cada lado. La referencia es el percentil 90 de todos los radios del vaso, y el estrechamiento es

```
s = min(1, max(0, 1 - r / r_ref))
```

Si `r` iguala a la referencia, `s` vale 0 (morado). Cuanto más angosto el vaso, más sube `s` hacia 1 (amarillo). La escala de color es fija entre 0 y 1, así que los casos son comparables. El tubo de 0,3 mm es solo visual y no representa el calibre.

## Resultados

| | `dataset08` | `dataset00` |
|---|---|---|
| Intensidad en la arteria | 345 | 324 |
| Triángulos (cruda, recorte, conectividad) | 595 348 / 46 923 / 38 058 | 556 569 / 28 499 / 18 109 |
| Radio mínimo y máximo (mm) | 0,548 y 3,094 | 0,259 y 2,908 |
| Referencia, percentil 90 (mm) | 1,470 | 1,128 |
| Mayor estrechamiento estimado | 0,627 | 0,771 |
| Punto de la línea central | 112 de 195 | 170 de 221 |

En `dataset00` el radio mínimo está por debajo del tamaño del vóxel (del orden de 0,4 mm), así que ese valor está en el límite de lo que el dato resuelve.

![Tubo coloreado en dataset08](https://github.com/user-attachments/assets/5e8b6727-38b0-48b2-89ce-b1608d4cee31)

## Referencias

- H. A. Kirişli et al. *Standardized evaluation framework for evaluating coronary artery stenosis detection, stenosis quantification and lumen segmentation algorithms in computed tomography angiography.* Medical Image Analysis 17(8), 859-876, 2013. <https://doi.org/10.1016/j.media.2013.05.007>
- M. Borkin et al. *Evaluation of artery visualizations for heart disease diagnosis.* IEEE TVCG 17(12), 2479-2488, 2011.
- C. A. Karlo et al. *A systematic approach for analysis, interpretation, and reporting of coronary CTA studies.* Insights into Imaging 3, 215-228, 2012.

## Autora

Sofia Guerra Jiménez. Computación Gráfica.

Los datos pertenecen a sus autores y no se redistribuyen aquí.
