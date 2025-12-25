# Trabajo_AS
Repositorio para el trabajo de la asignatura análisis de señales para el grupo C: Adrián Carrasco (Matemáticas), Clara Montalvá (ADE), Javier Herrero (Física), Fabián Calvo (Ingeniería Industrial), Juan Alcaraz (Matemáticas)

# Compresión de vídeo

## Compresión de imagen

## Compresión de audio
Esta sección se centra en la compresión de la señal sonora del vídeo utilizando la Transformada Wavelet Discreta (DWT). Para el desarrollo de este módulo se ha creado el script **`funciones_audio.R`**, disponible en la raíz de este repositorio.

### Metodología
La compresión se realiza en el dominio wavelet mediante el **umbralado selectivo de la energía** de los coeficientes. Los pasos a seguir son los siguientes:

1. **Descomposición Multirresolución**: Se aplican 6 niveles de descomposición, cubriendo las frecuencias medias y bajas donde se concentra la mayor densidad energética del audio (44100 Hz).
2. **Umbralado por Energía Acumulada**: Se ordenan los coeficientes por magnitud y se eliminan aquellos cuya contribución energética es inferior al umbral $\lambda$.
3. **Análisis de Entropía y Redundancia**: Se utiliza la **Entropía de Shannon** calculada a partir de histogramas de energía para medir la compresibilidad de la señal.

### Resultados
* **Punto Óptimo**: Se determinó que **$\lambda = 0.1$** es el valor de equilibrio ideal, logrando una alta reducción de memoria con una mínima pérdida de información relevante.
* **Eficiencia Energética**: Una fracción mínima de los coeficientes es capaz de capturar el 90% de la energía total del audio.
* **Formato de Salida**: Se utiliza el formato `.rds` para preservar la estructura de los coeficientes nulos, optimizando el almacenamiento real frente al formato `.wav`.
