Instrucciones del repositorio:

	El código para la guía está en 1 o más markdowns que hay que knitear a WORD (porque a odt queda feo el estilo del código).
	El word de cada sección (si fueran mas de 1) se incrusta en el master document 'Guia de TP RNAseq.odm', ahí está toda la info de los estilos del texto.
	Para esto: en el navegador del master document click derecho -> insertar -> buscas el subdocument (word en este caso)


Estructura del TP y comentarios:

	Introducción de anotación de genes
		ejemplos prácticos en el aula (con el web apollo abierto y metiendo mano)
		para ver intrones, exones, UTRs, splicing alternativo, etc...
			(buscar de antemano posiciones en las que hayan ejemplos concretos de cada cosa)
		BLAST y NCBI
		FASTAs
			
	Introducción al abordaje cuantitativo del análisis de RNAseq (mini-teórica que ya dí, mejorarla)
		tablas de conteos
		correcciones por longitud de c/gen y por reads totales  de la corrida # VOLARLO, no se implementa en DESeq
		comparación estadística para GENES particulares # buscar genes por id en i5k o apollo
		PCA, heatmaps, custering y similitud de los patrones de expresión entre MUESTRAS
		
		
	Cómo se hacen estas cosas en R
		Intro básica de elementos de R
			Objetos, funciones, argumentos, paquetes
		Guía (hay que curarla para sacarle complejidad, hay varias cosas innecesarias)
			Corregir strings de los filenames de entrada para matchear 'Conteos.csv' y 'Muestras.csv'
			


Planteamiento para hablar con Lore y Eze:
	Hacer primero un TP más clásico: 
		Una guia con pasos de un pipeline de ejemplo, se hace en al aula y se van viendo los problemas sobre la marcha (de R o conceptuales).
		En la guía de ejemplo se deberían mostrar:
			Comparación de expresión de genes puntuales entre contextos celulares
				Búsqueda de genes usando blast
				Identificación de esos genes en el dataset
				t-test (lo ven en estadística 1), gráficos custom, etc...		# OJO, ESTO NO ES TÉCNICAMENTE RIGUROSO PORQUE LOS DATOS VIENEN SIN TAMA NO DE LAS LIBRERÍAS RNAseq		
			Comparación de patrón de transcripción gral entre dos contextos celulares (tejidos, sexo y/o dsxRNai).
				PCA, heatmaps de muestras
		Esto se entregaría de forma sencilla. Se puede hacer un "ïnforma" que sea una guía de preguntas simples que irían contestando
			El objetivo del informe sería mas que nada para confirmar que estuvieron siguiendo el tp, nada más.
	Que después ellos (ellas en 2022) reproduzcan el pipeline que vimos juntos pero aplicado a una pregunta de su interés
		Este sería el TP final que presentan
		

		