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
├── data/                   # Carpeta donde se guardan las imágenes y los archivos de video
│   ├── frames/             # Carpeta donde se guardan las imágenes del vídeo (fotogramas)
│   ├── acantilado.jpg
│   ├── audio.wav
│   ├── barcos.png
│   ├── dibujo.png     
│   ├── flores.jpg     
│   ├── noche.jpg
│   ├── video.mp4
│   └── video_R.mp4    # Vídeo cargado, separado y directamente guardado
├── data_comp/                 # Resultados de la compresión
│   ├── frames_comp/         # Carpeta con las imágenes comprimidas del vídeo
│   ├── acantilado_DCT.jpg       # Imagen comprimida con DCT
│   ├── acantilado_DWT.jpg       # Imagen comprimida con DWT
│   ├── audio2comp.wav       # Audio comprimido
│   ├── audio_wavelet_comp.rds       # Archivo optimizado con ceros
│   ├── barcos_DCT.png       # Imagen comprimida con DCT
│   ├── barcos_DWT.png       # Imagen comprimida con DWT
│   ├── dibujo_DCT.png       # Imagen comprimida con DCT
│   ├── dibujo_DWT.png       # Imagen comprimida con DWT
│   ├── flores_DCT.jpg       # Imagen comprimida con DCT
│   ├── flores_DWT.jpg       # Imagen comprimida con DWT
│   ├── noche_DCT.jpg       # Imagen comprimida con DCT
│   ├── noche_DWT.jpg       # Imagen comprimida con DWT
│   └── video_comp.mp4         # Vídeo comprimido
├── explicaciones_y_pruebas/          # Scripts de código fuente
|   ├── IMAGENES SPIHT/      # Carpeta con imágenes comprimidas con SPIHT
|   ├── IMAGENES/            # Carpeta con imágenes de prueba para los algoritmos
|   ├── Captura de pantalla 2025-12-09 194338.png
│   ├── Explicacion_Algoritmos_DCT_IDWT.Rmd
│   ├── Explicacion_Algoritmos_EZW_SPIHT.Rmd
│   ├── Imagen_descomposicion_DWT.png
│   ├── Pruebas_Adrian.Rmd
│   ├── Pruebas_SPIHT.Rmd
│   ├── Reconstruccion_500_iter.jpeg
│   ├── Reconstruccion_50_iter.jpeg
│   └── audio_video.Rmd
├── funciones/      # Scripts de código fuente
│   ├── Codificacion_SPIHT_funcional.R   # Funciones para la codificación SPIHT
│   ├── Decodificacion_SPIHT.R
│   ├── compresion_DCT.R
│   ├── compresion_DCT_opt.R    # Función para compresión DCT
|   ├── compresion_DWT.R    # Función para compresión DWT
│   ├── compresion_IDWT.R
│   ├── compresion_SPIHT.R    # Funciones para la compresión SPIHT
│   ├── funciones_audio.R   # Funciones para el vídeo
│   └── metricas.R    # Funciones para las diferentes métricas
├── README.md      # Este archivo
├── Trabajo_AS.Rmd      # Script principal que ejecuta la compresión
├── Trabajo_AS.pdf      # Memoria
├── ieee.csl      # Memoria
└── Referencias.bib      # Referencias
```


## Módulos del proyecto

### Compresión DWT

### Compresión DCT

### Compresión EZW y SPIHT

Se explica de forma breve el funcionamiento de los algoritmos y se hace un ejemplo de compresión con una imagen de pequeñas dimensiones.

Después se comentan los resultados obtenidos y posibles causas de la falta de calidad, y posibles mejoras.

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
* **K.R. Rao & P.C. Yip (2001).** *The Transform and Data Compression Handbook* (1th ed.). CRC Press. (Capítulo 6: Wavelet Based Image Compression)
* **Sullivan, G. J., & Wiegand, T. (2005).** *Video compression—from concepts to the H.264/AVC standard.* *Proceedings of the IEEE*, 93(1), 18–31.
* **Perales Gómez, Á. L.** *DCT para la compresión de imágenes con perdida de calidad*. Documento Técnico. Universidad de Murcia, Facultad de Informática.
* **Wang, Z., Bovik, A. C., Sheikh, H. R., & Simoncelli, E. P. (2004).** *Image quality assessment: From error visibility to structural similarity.* *IEEE Transactions on Image Processing*, 13(4), 600–612.
