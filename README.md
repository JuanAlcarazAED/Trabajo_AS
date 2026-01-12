# Trabajo Análisis de Señales: Compresión de imagen y vídeo
Repositorio para el trabajo de la asignatura análisis de señales para el grupo C: 

Miembros del grupo: 
- Adrián Carrasco (Matemáticas)
- Clara Montalvá (ADE)
- Javier Herrero (Física)
- Fabián Calvo (Ingeniería Industrial)
- Juan Alcaraz (Matemáticas)


## Estructura del repositorio
```
Trabajo_AS/
├── im_2/                  # Carpeta donde se guardan las imágenes
│   ├── acantilado.jpg
│   ├── dibujo.jpg        
│   ├── flores.jpg     
│   └── noche.jpg
├── data/                   # Archivos de video
│   ├── video.mp4
│   ├── frames/             # Carpeta donde se guardan las imágenes
│   └── audio.wav
├── src/                    # Scripts de código fuente
│   ├── compresion_IDWT.R
│   ├── compresion_DCT.R
│   ├── compresion_DCT_opt.R
│   ├── metricas.R
│   ├── funciones_audio.R   # Funciones para el vídeo
│   ├── Codificacion_SPIHT.R
│   ├── compresion_SPIHT.R
│   └── Trabajo_AS.Rmd      # Script principal que ejecuta la compresión
├── data_comp/                 # Resultados de la compresión
│   ├── audio_wavelet_comp.rds       # Archivo optimizado con ceros 
│   ├── audio2comp.wav       # Audio comprimido
│   ├── audio_wavelet_comp.rds       # Archivo optimizado con ceros 
│   ├── frames_comp/         # Carpeta con las imagenes comprimidas del vídeo
│   ├── video_comp.mp4         # Vídeo comprimido
│   ├── acantilado_DCT.jpg       # Imagen comprimida con DCT
│   ├── dibujo_DCT.jpg       # Imagen comprimida con DCT
│   ├── flores_DCT.jpg       # Imagen comprimida con DCT
│   ├── noche_DCT.jpg       # Imagen comprimida con DCT
│   ├── acantilado_IDWT.jpg     # Imagen comprimida con IDWT
│   ├── dibujo_IDWT.jpg       # Imagen comprimida con IDWT
│   ├── flores_IDWT.jpg       # Imagen comprimida conIDWT
│   ├── noche_IDWT.jpg       # Imagen comprimida con IDWT
├── Trabajo_AS.pdf               # Memoria
├── Referencias.bib               # Referencias
```



## Módulos del proyecto

### Compresión IDWT

### Compresión DCT

### Compresión EZW y SPIHT


### Compresión de vídeo

Para la compresión del vídeo se va a comprimir la imagen del vídeo y el audio por separado.

#### Compresión de imagen

#### Compresión de audio
Esta sección se centra en la compresión de la señal sonora del vídeo utilizando la Transformada Wavelet Discreta (DWT). Para el desarrollo de este módulo se ha creado el script **`funciones_audio.R`**, disponible en la raíz de este repositorio.

La compresión se realiza en el dominio wavelet mediante el **umbralado selectivo de la energía** de los coeficientes. Los pasos a seguir son los siguientes:

1. **Descomposición Multirresolución**: Se aplican 6 niveles de descomposición, cubriendo las frecuencias medias y bajas donde se concentra la mayor densidad energética del audio (44100 Hz).
2. **Umbralado por Energía Acumulada**: Se ordenan los coeficientes por magnitud y se eliminan aquellos cuya contribución energética es inferior al umbral $\lambda$.
3. **Análisis de Entropía y Redundancia**: Se utiliza la **Entropía de Shannon** calculada a partir de histogramas de energía para medir la compresibilidad de la señal.

En cuanto a los resultados obtenidos

* **Punto Óptimo**: Se determinó que **$\lambda = 0.1$** es el valor de equilibrio ideal, logrando una alta reducción de memoria con una mínima pérdida de información relevante.
* **Eficiencia Energética**: Una fracción mínima de los coeficientes es capaz de capturar el 90% de la energía total del audio.
* **Formato de Salida**: Se utiliza el formato `.rds` para preservar la estructura de los coeficientes nulos, optimizando el almacenamiento real frente al formato `.wav`.

## Instrucciones de uso

Para reproducir el trabajo, ejecute el archivo `Trabajo_AS.Rmd`

## Referencias

* **Mallat, S. (2009).** *A Wavelet Tour of Signal Processing: The Sparse Way* (3rd ed.). Academic Press. (Capítulos 1–6).
* **Gonzalez, R. C., & Woods, R. E. (2018).** *Digital Image Processing* (4th ed.). Pearson. (Capítulo 8: Image Compression).
* **Sullivan, G. J., & Wiegand, T. (2005).** Video compression—from concepts to the H.264/AVC standard. *Proceedings of the IEEE*, 93(1), 18–31.
* **Perales Gómez, Á. L.** *DCT para la compresión de imágenes con perdida de calidad*. Documento Técnico. Universidad de Murcia, Facultad de Informática.

