"""_summary_
Dado el siguiente mensaje, emitido por una fuente de información:

Obtener las probabilidades de los símbolos y la entropía de la fuente
Generar un código binario óptimo para la fuente de información
Determinar la tasa de compresión que se obtiene codificando el mensaje
Generar un código binario óptimo para la extensión de orden 3
Calcular la longitud media, el rendimiento y la redundancia de cada código
Verificar si ambos códigos cumplen con el Primer Teorema de Shannon
"""
import math

def alfabeto(texto: str):
    """_summary_
    Dada una cadena de caracteres, devuelve una lista con su alfabeto y una lista con sus
    probabilidades.

    Args:
        texto (str): La cadena de caracteres.
    Returns:
        alf (list), probabilidades (list): alfabeto de la cadena, probabilidades de cada símbolo
        del alfabeto.
    """
    longitud = len(texto)
    alf = list(set(texto))  # alfabeto como lista desde un set
    probabilidades = [texto.count(simbolo) / longitud for simbolo in alf]

    return alf, probabilidades

def alfabeto_codigo(codigo:list) -> str:
    """_summary_
        Dado un código, devuelve una cadena con el alfabeto del código.
    Args:
        codigo (list): Código
    Returns:
        str: Alfabeto del código
    """
    alfabeto = set()
    for palabra in codigo:
        for char in palabra:
            alfabeto.add(char)
    return ''.join(sorted(alfabeto))

def cantidad_informacion(p: float, r=2) -> float:
    """_summary_
    Calcula la cantidad de información de un parto dado su probabilidad.

    Args:
        p (float): La probabilidad del parto.
        r (int): Unidad de medida de la informacion.
    Returns:
        resultado (float): La cantidad de información en bits.
    """
    try:
        if p > 1 or p <= 0:
            raise ValueError("Probabilidad no válida")
        resultado = math.log(1/p, r)
    except ValueError:
        resultado = -1.0
    return resultado

def entropia(probabilidades: list, r=2) -> float:
    """_summary_
    Calcula la entropía de una fuente de información dada una lista de probabilidades.

    Args:
        probabilidades (list): Lista de probabilidades de los partos.
        r (int): Unidad de medida de la informacion.
    Returns:
        ent (float): La entropía en la base especificada
    """
    ent = 0
    for p in probabilidades:
        ent += p * cantidad_informacion(p, r)
    return ent

def entropia_codigo(codigo:list, probabilidades:list) -> float:
    """_summary_
        Calcula la entropía de un código dado sus probabilidades.
    Args:
        codigo (list): Código
        probabilidades (list): Probabilidades de cada símbolo en el código
    Returns:
        float: Entropía del código
    """
    num_simbolos = len(alfabeto_codigo(codigo))

    return entropia(probabilidades, r=num_simbolos)

def longitudes(codigo:list) -> list:
    """_summary_
        Dado un código, devuelve una lista con las longitudes de las palabras del código.
    Args:
        codigo (list): Código
    Returns:
        list: Longitudes de las palabras del código
    """
    return [len(palabra) for palabra in codigo]

def longitud_media(codigo:list, probabilidades:list) -> float:
    """_summary_
        Calcula la longitud media de un código dado sus probabilidades.
    Args:
        codigo (list): Código
        probabilidades (list): Probabilidades de cada símbolo en el código
    Returns:
        float: Longitud media del código
    """
    long = longitudes(codigo)
    return sum(l * p for l, p in zip(long, probabilidades))

def generar_combinaciones(alfabeto:list, orden:int):
    """_summary_
        Genera todas las combinaciones posibles de un alfabeto dado un orden
        Recursivo
    Args:
        alfabeto (list): Alfabeto de la fuente
        orden (int): Orden de la extensión
    Returns:
        dict: Diccionario con las combinaciones posibles, y su composición
        ejemplo: {"BABABA": ["BA","BA","BA"]...}. Sirve para cualquier longitud de palabra :)
    """
    if orden == 1:
        combinaciones = {}
        for simbolo in alfabeto:
            combinaciones[simbolo] = [simbolo]
    else:
        combinaciones = {}
        for simbolo in alfabeto:
            combinaciones_ant = generar_combinaciones(alfabeto, orden - 1)
            for key,value in combinaciones_ant.items():
                combinaciones[simbolo + key] = [simbolo]
                for combinacion in value:
                    combinaciones[simbolo + key].append(combinacion)
    return combinaciones

def extension_orden_n(alfabeto:list, probabilidades:list, orden:int):
    """_summary_
        Devuelve la extension de orden N de una fuente
    Args:
        alfabeto (list): Simbolos de la fuente
        probabilidades (list): Probabilidades de cada simbolo
        n (int): Orden de la extension
    """
    extension = generar_combinaciones(alfabeto, orden)
    probs = []
    for key,value in extension.items():
        prob_palabra = 1
        for char in value:
            prob_palabra *= probabilidades[alfabeto.index(char)]
        probs.append(prob_palabra)
    return extension, probs

def verifica_shannon(probabilidades:list, palabras_codigo:list, n:int):
    """_summary_
        Verifica el teorema de Shannon para la extension de orden n de una fuente de información
    Args:
        probabilidades (list): Lista de probabilidades de los símbolos de la fuente
        palabras_codigo (list): Palabras codigo
        n (int): Orden de la extensión
    """
    if n > 1:
        ext,probs = extension_orden_n(palabras_codigo, probabilidades, n)
        ext = list(ext.keys())
    else:
        ext = palabras_codigo
        probs = probabilidades
    long_med = longitud_media(ext, probs)
    ent = entropia_codigo(palabras_codigo, probabilidades)
    print(f"Entropía: {ent:.4f}, Longitud media/n: {long_med/n:.4f}")
    cumple = ent <= long_med/n <= ent + 1/n
    return cumple


def huffman_recursivo(prob_dict: dict) -> dict:
    """
    Función recursiva auxiliar para huffman
    """
    # CASO BASE, devuelvo 0 y 1 para los dos símbolos
    if len(prob_dict) == 2:
        simbolos = list(prob_dict.keys())
        return {simbolos[0]: '0', simbolos[1]: '1'}

    # RECURSIVIDAD
    else:
        diccionario_ordenado = dict(sorted(prob_dict.items(), key=lambda item: item[1]))

        # Combino los simbolos de menor probabilidad y armo un alfabeto reducido
        simbolos = list(diccionario_ordenado.keys())
        simbolo1 = simbolos[0]
        simbolo2 = simbolos[1]

        prob_simbolo1 = diccionario_ordenado[simbolo1]
        prob_simbolo2 = diccionario_ordenado[simbolo2]

        nuevo_simbolo = simbolo1 + simbolo2
        nueva_probabilidad = prob_simbolo1 + prob_simbolo2

        diccionario_reducido = {k: v for k, v in diccionario_ordenado.items() if k not in [simbolo1, simbolo2]}
        diccionario_reducido[nuevo_simbolo] = nueva_probabilidad

        codigo_reducido = huffman_recursivo(diccionario_reducido)

        # Obtengo el código asignado al nuevo símbolo
        codigo_para_nuevo_simbolo = codigo_reducido[nuevo_simbolo]

        # Asigno el valor a los cods originales y los diferencio con un 0 y un 1
        codigo_reducido[simbolo1] = codigo_para_nuevo_simbolo + "0"
        codigo_reducido[simbolo2] = codigo_para_nuevo_simbolo + "1"

        # Elimino el simbolo combinado
        del codigo_reducido[nuevo_simbolo]

        return codigo_reducido

def huffman(probabilidades: list) -> dict:
    """
        Calcula un código Huffman óptimo dada una lista de probabilidades.
    Args:
        probabilidades (list): Lista de probabilidades de los símbolos.
    Returns:
        dict: Código óptimo. Las claves son los índices de la lista
              original y los valores son los códigos binarios.
              Ej: {0: '0', 1: '10', 2: '11'}
    """
    if not probabilidades:
        return {}

    if len(probabilidades) == 1:
        return {0: '0'}

    # Lista de probs en forma de diccionario
    prob_dict = {str(i): p for i, p in enumerate(probabilidades)}

    codigo_optimo = huffman_recursivo(prob_dict)

    # Convertir las claves de str a int (No es necesario pero es más limpio)
    codigo_final = {int(k): v for k, v in codigo_optimo.items()}

    return codigo_final

def sf_recursivo(probabilidades: list, codigo_actual: str, codigos: dict):
    """
    Función recursiva auxiliar para shannon_fano
    """

    if len(probabilidades) == 1:
        original_index = probabilidades[0][1]
        codigos[original_index] = codigo_actual or "0"
        return

    total = sum(p for p, _ in probabilidades)
    mitad = total / 2.0

    acum = 0
    corte = 0

    for i, (prob, _) in enumerate(probabilidades):
        acum += prob

        # Si alcanzamos o superamos la mitad, evaluamos si cortar aquí o antes
        if acum >= mitad:

            # Calcular la diferencia si incluimos el elemento actual en el Grupo 1
            diff_actual = abs(acum - mitad)

            # Calcular la diferencia si no incluimos el actual
            suma_prev = acum - prob
            diff_sin_actual = abs(suma_prev - mitad)

            # Si la suma anterior estaba más cerca de la mitad el indice de corte es i-1
            if diff_sin_actual < diff_actual:
                corte = i
            else:
                corte = i + 1
            break

    if corte < 1:
        corte = 1
    elif corte >= len(probabilidades):
        corte = len(probabilidades) - 1

    group1 = probabilidades[:corte]
    group2 = probabilidades[corte:]

    sf_recursivo(group1, codigo_actual + "0", codigos)
    sf_recursivo(group2, codigo_actual + "1", codigos)


def shannon_fano(probabilidades: list) -> dict:
    """
        Calcula un código Shannon-Fano dada una lista de probabilidades.
    Args:
        probabilidades (list): Lista de probabilidades de los símbolos.
                               Ej: [0.5, 0.2, 0.3]
    Returns:
        dict: Código. Las claves son los índices de la lista
              original y los valores son los códigos binarios.
              Ej: {0: '0', 1: '10', 2: '11'}
    """
    # Indexo las probabilidades y ordeno descendente
    probs = [(p, i) for i, p in enumerate(probabilidades)]
    probs.sort(key=lambda x: x[0], reverse=True)

    codigos = {}

    if probs:
        sf_recursivo(probs, "", codigos)

    return codigos


def rendimiento_y_redundancia(probabilidades: list, codigo: list) -> tuple[float, float]:
    """
    Calcula el rendimiento y la redundancia de un código dado sus probabilidades.

    Args:
        probabilidades (list): Lista de probabilidades de los símbolos.
        codigo (list): Lista de palabras código.
    Returns:
        tuple: Rendimiento y redundancia del código.
    """
    ent = entropia_codigo(codigo, probabilidades)
    long_med = longitud_media(codigo, probabilidades)
    print(ent,long_med)
    rendimiento = ent / long_med if long_med != 0 else 0
    return rendimiento, 1 - rendimiento

def rendimiento_redundancia_2(ent, long_media):
    """_summary_
        Reversion REND Y REDU
    Args:
        entropia (_type_): _description_
        longitud_media (_type_): _description_
    """
    if long_media != 0:
        rendimiento = ent/long_media
    else:
        rendimiento = 0
    return rendimiento, 1 - rendimiento


def codificar(mensaje: str, alfabeto: list[str], codificacion: list[str]) -> bytearray:
    """
    Codifica un mensaje string a un bytearray usando una tabla de codigos.
    El primer byte contiene en sus 3 bits más significativos la cantidad de bits
    de relleno (residuo) que hay en el último byte.
    """
    if not mensaje:
        raise ValueError("El mensaje no puede estar vacío.")

    mensaje_codificado = bytearray()
    #Uso los 3 primeros bits para contabilizar el residuo, inicializo en 0
    buffer_bits = "000"

    for simbolo in mensaje:
        if simbolo in alfabeto:
            indice = alfabeto.index(simbolo)
            buffer_bits += codificacion[indice]

            # Cada vez que juntamos 8 bits, creamos un byte
            while len(buffer_bits) >= 8:
                byte_val = int(buffer_bits[:8], 2)
                mensaje_codificado.append(byte_val)
                buffer_bits = buffer_bits[8:]
        else:
            raise ValueError(f"El símbolo '{simbolo}' no está en el alfabeto fuente.")

    # Calculo que residuo
    residuo = 8 - len(buffer_bits) if buffer_bits else 0

    # Rellenamos con ceros a la derecha para completar el octeto
    mensaje_codificado.append(int(buffer_bits.ljust(8, '0'), 2))

    # Ingreso el residuo en los 3 bits más significativos del primer byte
    mensaje_codificado[0] = (residuo << 5) | mensaje_codificado[0]

    return mensaje_codificado

def tasa_de_compresion(mensaje_original: str, mensaje_codificado: bytearray) -> float:
    """
    Calcula la tasa de compresión entre el mensaje original y el mensaje codificado.
    Args:
        mensaje_original (str): El mensaje original en formato string.
        mensaje_codificado (bytearray): El mensaje codificado en formato bytearray.
    Returns:
        float: La tasa de compresión.
    """
    tamano_original = len(mensaje_original) * 8  # en bits
    tamano_codificado = len(mensaje_codificado) * 8  # en bits

    if tamano_codificado == 0:
        raise ValueError("El mensaje codificado no puede tener tamaño cero.")

    tasa_compresion = tamano_original / tamano_codificado

    return tasa_compresion

MENSAJE = "GIIHGGGHGIHHIHIIGFHH"

if __name__ == "__main__":

    alfabeto_msj,probs = alfabeto(MENSAJE)
    ent = entropia(probs)
    print(f"Entropia de la fuente: {ent:.4f}")
    dict_huffman = huffman(probs)
    codigo_huffman = [dict_huffman[i] for i in range(len(probs))]
    dict_sf = shannon_fano(probs)
    codigo_sf = [dict_sf[i] for i in range(len(probs))]

    rend_huffman, redu_huffman = rendimiento_y_redundancia(probs,codigo_huffman)
    rend_sf, redu_sf = rendimiento_y_redundancia(probs,codigo_sf)

    if rend_huffman > rend_sf:
        codigo_fuente = codigo_huffman
        codificacion = codificar(MENSAJE,alfabeto_msj,codigo_huffman)
        rendimiento_fuente = rend_huffman
        redundancia_fuente = redu_huffman
    else:
        codigo_fuente = codigo_sf
        codificacion = codificar(MENSAJE,alfabeto_msj,codigo_sf)
        rendimiento_fuente = rend_sf
        redundancia_fuente = redu_sf
    print(f"Codigo fuente: {codigo_fuente}")
    tdc = tasa_de_compresion(MENSAJE, codificacion)
    print(f"Tasa de compresion de la codificacion: {tdc:.4f}")
    ####### EXTENSION ORDEN 3

    ext_ord3, probs_ord3 = extension_orden_n(alfabeto_msj,probs,3)
    dict3_huffman = huffman(probs_ord3)
    ord3_huffman = dict3_huffman.values()
    dict3_sf = shannon_fano(probs_ord3)
    ord3_sf = [dict3_sf[i] for i in range(len(probs_ord3))]

    rend3_huffman, redu3_huffman = rendimiento_y_redundancia(probs_ord3,codigo_huffman)
    rend3_sf, redu3_sf = rendimiento_y_redundancia(probs_ord3,codigo_sf)

    if rend3_huffman > rend3_sf:
        codigo_ord3 = ord3_huffman
        rendimiento_ord3 = rend3_huffman
        redundancia_ord3 = redu3_huffman
    else:
        codigo_ord3 = ord3_sf
        rendimiento_ord3 = rend3_sf
        redundancia_ord3 = redu3_sf
    print(f"Codigo de la extension orden 3: {codigo_ord3}")
    ##########

    long_media_fuente = longitud_media(codigo_fuente, probs)
    long_media_ord3 = longitud_media(codigo_ord3, probs_ord3)

    rend3v2, redu3v2 = rendimiento_redundancia_2(entropia_codigo(codigo_ord3,probs_ord3),long_media_ord3)

    shannon_fuente = verifica_shannon(probs, codigo_fuente,1)
    shannon_ord3 = verifica_shannon(probs_ord3, codigo_ord3,3)

    print("Fuente:\n")
    print(f"Longitud media: {long_media_fuente:.4f}")
    print(f"Rendimiento: {rendimiento_fuente:.4f}")
    print(f"Redundancia: {redundancia_fuente:.4f}")
    print(f"Verifica Shannon? {shannon_fuente}")

    print()

    print("Extension de orden 3\n")
    print(f"Longitud media: {long_media_ord3:.4f}")
    print(f"Rendimiento: {rend3v2:.4f}")
    print(f"Redundancia: {redu3v2:.4f}")
    print(f"Verifica Shannon? {shannon_ord3}")



