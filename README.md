Instrucciones del repositorio:

	El código para la guía está en 1 o más markdowns que hay que knitear a WORD (porque a odt queda feo el estilo del código).
	El word de cada sección (si fueran mas de 1) se incrusta en el master document 'Guia de TP RNAseq.odm', ahí está toda la info de los estilos del texto.
	Para esto: en el navegador del master document click derecho -> insertar -> buscas el subdocument (word en este caso)


Estructura del TP y comentarios:

	Introducción de anotación de genes
		Repaso de RNA-seq -> Presentación Licán.
		Ver intrones, exones, UTRs, splicing alternativo, etc... (buscar genes piolas si hay tiempo)
		BLAST y NCBI
		FASTAs 
			
	Introducción al abordaje cuantitativo del análisis de RNAseq (mini-teórica que ya dí, mejorarla)
	(seguir presentación de lican)
		tablas de conteos
		comparación estadística para GENES particulares # buscar genes por id en i5k o apollo
		PCA, heatmaps, custering y similitud de los patrones de expresión entre MUESTRAS
		
		
	Cómo se hacen estas cosas en R
		Intro básica de elementos de R (ANEXO)
			Objetos, funciones, argumentos, paquetes (ANEXO)
	Guía Principal (hay que curarla para sacarle complejidad, hay varias cosas innecesarias)
		
			


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
		


Cronograma final de actividades:

	1. 6 DE MAYO: TP4a (PRE-TP4a) -> Transcriptómica y anotación (2 hs) (OFFLINE, EN AULA 9)
		1.1 Introducción a transcriptomica:
			- Teórica RNA-sec repaso ilúmina.
			- Visualización para ver intrones, exones, UTRs, splicing alternativo, etc... (buscar de antemano posiciones en las que hayan ejemplos concretos de cada cosa EN APOLLO)
		1.2 BLAST y NCBI
		1.3 FASTAs 
		1.4 Introducción al abordaje cuantitativo del análisis de RNAseq (mini-teórica que ya dí, mejorarla)
			- tablas de conteos
			- correcciones por longitud de c/gen y por reads totales  de la corrida # VOLARLO, no se implementa en DESeq
			- comparación estadística para GENES particulares # buscar genes por id en i5k o apollo
			- PCA, heatmaps, custering y similitud de los patrones de expresión entre MUESTRAS


	2. 13 DE MAYO: TP4b (TP4-1° parte) -> Genómica II (POCO INTERNET, EN AULA 9 o enroque con el miércoles anterior.) TENER LA GUÍA LISTA PARA ESTA FECHA.
		2.1 Teórica edu sobre los ontophagus taurus ? Si nó que la miren grabada. (suma el RNAi de dsx) (1 h + 15' de recreo para recuperarse del impacto)
		2.2 Cómo se hacen estas cosas en R (ANEXO) (últimos 45' minutos)
			Intro básica de elementos de R (ANEXO)
			Objetos, funciones, argumentos, paquetes (ANEXO)
		
			

	3. 27 DE MAYO: TP5 (TP4-2° parte)-> RNA-sec II + planteo del TP final. (MUCHO INTERNET, AULA A CONFIRMAR-POSÍBLEMENTE EN LABO BIO I O II)
		3.1 Arrancar Guía análisis de RNA-sec(hay que curarla para sacarle complejidad, hay varias cosas innecesarias) - que vayan respondiendo preguntas presentes en la Guía.
			3.1.1 VISUALIZACION
			3.1.2 DE analisis clásico
				3.1.2.1 Análisis de DE 
				3.1.2.2 buscar genes en NCBI. (Posiblemente dárselos de tarea)
			3.1.3 ir a buscar un gen conocido para ver si tiene DE ??
		3.2 Proponer trabajo final.


Materiales para las chicas:

	1. Teórica de Lican.
	3. Guía análisis de RNA-sec.
	2. Guía "complementaria" ? "Guía R + conceptos básicos en bioinformática: BLAST, GEO, NCBI, FASTAs, etc".