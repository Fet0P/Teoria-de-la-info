"""_summary_
Dadas las siguientes secuencias de entrada y salida de un canal:

Obtener las probabilidades a priori y la matriz del canal
Informar si se trata de un canal sin ruido y/o determinante
Calcular la equivocación, la pérdida y la información mutua
Efectuar todas las reducciones suficientes posibles en el canal
Determinar la capacidad y la probabilidad de error del canal
Detectar y/o corregir los errores en la secuencia de entrada
Aclaración: asumir que la secuencia de entrada contiene un mensaje representado con código ASCII 
y sus bits de paridad vertical, longitudinal y cruzada.
"""

import math

def imprimir_matriz(matriz, filas, cols):
    """_summary_
        Imprime una matriz con sus filas y columnas etiquetadas.
    """
    header = "      " + "  ".join([f"{c:>5}" for c in cols])
    print(header)
    print("-" * len(header))

    for i, fila_sym in enumerate(filas):
        valores = "  ".join([f"{val:5.4f}" for val in matriz[i]])
        print(f"{fila_sym:>3} | {valores}")

def generar_matriz_canal(str_entrada, str_salida):
    """
    Genera una matriz de transición P(b|a)
    Args:
        str_entrada (str): La cadena enviada al canal.
        str_salida (str): La cadena recibida del canal.
    Retorna:
        matriz (list of lists): Los valores de probabilidad.
        filas (list): Los símbolos de entrada (a).
        columnas (list): Los símbolos de salida (b).
    """

    simbolos_x = sorted(list(set(str_entrada)))
    simbolos_y = sorted(list(set(str_salida)))

    n_filas = len(simbolos_x)
    n_cols = len(simbolos_y)

    map_a = {char: i for i, char in enumerate(simbolos_x)}
    map_b = {char: i for i, char in enumerate(simbolos_y)}

    matriz = [[0.0] * n_cols for _ in range(n_filas)]
    total_por_fila = [0] * n_filas

    for a, b in zip(str_entrada, str_salida):
        fila_idx = map_a[a]
        col_idx = map_b[b]

        matriz[fila_idx][col_idx] += 1
        total_por_fila[fila_idx] += 1

    for i in range(n_filas):
        total = total_por_fila[i]
        if total > 0:
            for j in range(n_cols):
                matriz[i][j] = matriz[i][j] / total

    return matriz, simbolos_x, simbolos_y

def calcular_a_priori(entrada:str):
    """
    Calcula las probabilidades a priori P(a) de los símbolos de entrada.
    Args:
        str_entrada (str): La cadena enviada al canal.
    Returns:
        dict: Un diccionario donde la clave es el símbolo y
              el valor es su probabilidad.
    """
    longitud_total = len(entrada)

    if longitud_total == 0:
        return {}

    conteo = {}
    for char in entrada:
        if char in conteo:
            conteo[char] += 1
        else:
            conteo[char] = 1

    probs = {}
    for char in sorted(conteo.keys()):
        probs[char] = conteo[char] / longitud_total

    return probs

def es_canal_sin_ruido(matriz_canal:list[list[float]]) -> bool:
    """_summary_
        Verifica si el canal es sin ruido
        En la matriz: Ninguna columna debe tener más de una entrada distinta de cero.
    Args:
        matriz_canal (list[list[float]]): Matriz de probabilidades condicionales P(B|A)
    Returns:
        bool: True si el canal es sin ruido, False en caso contrario
    """
    rows = len(matriz_canal)
    cols = len(matriz_canal[0])

    for j in range(cols):
        entradas_a_j = 0

        for i in range(rows):
            if matriz_canal[i][j] < 1e-9: #Por trabajar con floats
                entradas_a_j += 1

        if entradas_a_j > 1:
            return False

    return True

def es_canal_determinante(matriz_canal:list[list[float]]) -> bool:
    """_summary_
        Verifica si el canal es determinante
        En la matriz: Ninguna fila debe tener más de una entrada distinta de cero.
    Args:
        matriz_canal (list[list[float]]): Matriz de probabilidades condicionales P(B|A)
    Returns:
        bool: True si el canal es determinante, False en caso contrario
    """
    for fila in matriz_canal:
        tiene_uno = False
        for valor in fila:
            if abs(valor - 1.0) < 1e-9: #Por trabajar con floats
                tiene_uno = True
                break

        if not tiene_uno:
            return False

    return True

def calcular_ruido(probs_priori:list[float],matriz_canal:list[list[float]]):
    """
    Calcula el ruido H(A|B) por definición.

    Args:
        probs_priori (list): Lista de probabilidades a priori P(a).
        matriz_canal (list of lists): Matriz de transición P(B|A).

    Returns:
        float: El ruido en bits.
    """
    rows = len(matriz_canal)
    cols = len(matriz_canal[0])

    probs_b = [0.0] * cols
    for j in range(cols):
        suma_col = 0.0
        for i in range(rows):
            suma_col += probs_priori[i] * matriz_canal[i][j]
        probs_b[j] = suma_col

    ruido = 0.0

    for i in range(rows):
        for j in range(cols):

            # Calculamos P(a, b) = P(a) * P(b|a)
            prob_conjunta = probs_priori[i] * matriz_canal[i][j]

            if prob_conjunta > 0:

                # P(b|a) = P(a, b) / P(b)
                prob_posteriori = prob_conjunta / probs_b[j]

                # H(A|B) -= P(a,b) * log2( P(b|a) )
                if prob_posteriori > 0:
                    ruido -= prob_conjunta * math.log2(prob_posteriori)

    return ruido

def calcular_perdida(probs_priori:list[float],matriz_canal:list[list[float]]):
    """
    Calcula la pérdida H(B|A)
    Args:
        probs_priori (list): Lista de probabilidades a priori P(a).
        matriz_canal (list of lists): Matriz de transición P(B|A).

    Returns:
        float: La perdida en bits.
    """
    perdida_total = 0.0

    for i, fila in enumerate(matriz_canal):
        p_a = probs_priori[i]

        # Calculamos la entropía de esta fila específica: H(B | A=ai)
        entropia_fila = 0.0
        for p_cond in fila:
            if p_cond > 0:
                entropia_fila -= p_cond * math.log2(p_cond)

        # Ponderamos por la probabilidad de que ocurra esta entrada
        perdida_total += p_a * entropia_fila

    return perdida_total

def calcular_informacion_mutua(probs_priori:list[float],matriz_canal:list[list[float]]):
    """
    Calcula la Información Mutua I(A;B) del canal.
    Args:
        probs_priori (list): Probabilidades a priori P(A).
        matriz_canal (list of lists): Matriz de transición P(B|A).

    Returns:
        float: Información mutua en bits.
    """
    rows = len(matriz_canal)
    cols = len(matriz_canal[0])

    # Calcular P(B)
    probs_b = [0.0] * cols
    for j in range(cols):
        suma_b = 0.0
        for i in range(rows):
            suma_b += probs_priori[i] * matriz_canal[i][j]
        probs_b[j] = suma_b

    info_mutua = 0.0

    for i in range(rows):
        p_a = probs_priori[i]

        for j in range(cols):
            p_b_dado_a = matriz_canal[i][j] # P(b|a)
            p_b = probs_b[j]                # P(b)

            if p_a > 0 and p_b_dado_a > 0 and p_b > 0:

                # P(a,b) * log2( P(b|a) / P(b) )
                # con P(a,b) = p_a * p_b_dado_a
                arg_log = p_b_dado_a / p_b
                prob_conjunta = p_a * p_b_dado_a

                info_mutua += prob_conjunta * math.log2(arg_log)

    return info_mutua

def multiplicar_matrices(matA:list[list[float]], matB:list[list[float]]):
    """
    Multiplica dos matrices
    Args:
        matriz_1 (list[list[float]]): Matriz P(A|B) de tamaño [N x M].
        matriz_2 (list[list[float]]): Matriz P(B|C) de tamaño [M x P].
    Returns:
        list[list[float]]: Matriz resultante P(A|C) de tamaño [N x P].
    """
    rows_1 = len(matA)
    cols_1 = len(matA[0])

    rows_2 = len(matB)
    cols_2 = len(matB[0])

    if cols_1 != rows_2:
        raise ValueError("Dimensiones incompatibles")

    matriz_resultante = [[0.0] * cols_2 for _ in range(rows_1)]

    for i in range(rows_1):
        for k in range(cols_2):
            suma = 0.0
            for j in range(cols_1):
                val_1 = matA[i][j]
                val_2 = matB[j][k]
                suma += val_1 * val_2

            matriz_resultante[i][k] = suma

    return matriz_resultante

def verificar_reduccion(matriz:list[list[float]], c1:int, c2:int):
    """
    Verifica si dos columnas de una matriz de canal son proporcionales.
    Si lo son, combinarlas representa una 'reducción suficiente' (sin pérdida de información).

    Args:
        matriz (list of lists): La matriz P(B|A).
        c1 (int): Índice de la primera columna.
        c2 (int): Índice de la segunda columna.
    Returns:
        bool: True si se pueden combinar, False en caso contrario.
    """
    ratio = None
    for fila in matriz:
        v1, v2 = fila[c1], fila[c2]
        if abs(v2) < 1e-9:
            if abs(v1) > 1e-9:
                return False
            continue
        r = v1 / v2
        if ratio is None:
            ratio = r
        elif abs(r - ratio) > 1e-5:
            return False
    return True

def generar_matriz_determinante(matriz: list[list[float]], c1: int, c2: int):
    """
    Genera una matriz de determinante para la reducción de canales
    utilizando listas estándar de Python.
    Args:
        matriz (list[list[float]]): La matriz original.
        c1 (int): El índice de la columna DESTINO (donde se suma).
        c2 (int): El índice de la columna ORIGEN (la que se elimina).

    Returns:
        list[list[int]]: La matriz determinante
    """
    n_cols = len(matriz[0])
    n_nuevas = n_cols - 1
    T = [[1 if x==y else 0 for y in range(n_nuevas)] for x in range(n_nuevas)]
    fila_nueva = [0] * n_nuevas
    idx_dest = c1 if c1 < c2 else c1 - 1
    fila_nueva[idx_dest] = 1
    T.insert(c2, fila_nueva)
    return T

def reducir_matriz(matriz: list[list[float]]) -> list[list[float]]:
    """_summary_
        Reduce la matriz lo maximo posible. Multiplicando sucesivamente la reduccion 
        por la matriz determinante correspondiente.
    Args:
        matriz (list[list[float]]): La matriz original.
    Returns:
        list[list[f.oat]]: La matriz determinante
    """
    matriz_reducida = [fila[:] for fila in matriz]

    i = 0
    j = 1

    if not matriz_reducida:
        return []
    n = len(matriz_reducida[0])

    while i < n and j < n:

        if verificar_reduccion(matriz_reducida, i, j):
            matriz_determinante = generar_matriz_determinante(matriz_reducida, i, j)
            matriz_reducida = multiplicar_matrices(matriz_reducida, matriz_determinante)
            n -= 1
            i = 0
            j = 1
        else:
            j += 1
            if j == n:
                i += 1
                j = i + 1

    return matriz_reducida

def es_canal_uniforme(matriz_canal:list[list[float]]):
    """
    Un canal es uniforme si cada fila de su matriz es una permutación de la primera fila

    Args:
        matriz_canal (list of lists): Matriz P(B|A).
    Returns:
        bool: True si todas las filas tienen los mismos valores (en cualquier orden).
    """
    if not matriz_canal:
        return False

    firma_referencia = sorted(matriz_canal[0])

    for i in range(1, len(matriz_canal)):
        fila_actual_ordenada = sorted(matriz_canal[i])

        if len(fila_actual_ordenada) != len(firma_referencia):
            return False

        for val_ref, val_act in zip(firma_referencia, fila_actual_ordenada):
            if abs(val_ref - val_act) > 1e-9:
                return False

    return True

def capacidad_max_info(matriz_canal:list[list[float]], step:float):
    """_summary_
        Devuelve la capacidad de un canal binario junto a la probabilidad óptima de entrada.
    Args:
        matriz_canal (list[list[float]]): Matriz BINARIA de transición P(B|A).
        step (float): Paso
    Returns:
        tuple: Capacidad del canal y probabilidad óptima de entrada.
    """
    if len(matriz_canal) != 2:
        raise ValueError("El canal no es binario")

    capacidad = -1.0
    probabilidad_optima = 0.0
    delta = 0.0

    while delta <= 1.0 + 1e-9:

        if delta > 1.0:
            delta = 1.0

        p_apriori = [delta, 1.0 - delta]

        info_mutua = calcular_informacion_mutua(p_apriori, matriz_canal)

        if info_mutua > capacidad:
            capacidad = info_mutua
            probabilidad_optima = delta

        delta += step

    return probabilidad_optima, capacidad

def prob_error_canal(probs_a_priori: list[float], matriz_canal: list[list[float]]):
    """_summary_
        Calcula la probabilidad de error de un canal
    Args:
        probs_a_priori (list[float]): Probabilidades a priori de las entradas.
        matriz_canal (list[list[float]]): Matriz de transición P(B|A).

    Returns:
        float: La probabilidad de error del canal
    """
    rows = len(matriz_canal)
    cols = len(matriz_canal[0])

    regla_de_decision = []
    for j in range(cols):
        max_val = -1.0
        mejor_input = -1
        for i in range(rows):
            if matriz_canal[i][j] > max_val:
                max_val = matriz_canal[i][j]
                mejor_input = i
        regla_de_decision.append(mejor_input)

    prob_error = 0.0
    for i in range(rows):
        for j in range(cols):
            if regla_de_decision[j] != i:
                prob_error += probs_a_priori[i] * matriz_canal[i][j]

    return prob_error

def analizar_paridad(matriz: list[str]) -> dict:
    """_summary_
        Analiza una matriz de mensaje
        Devuelve un diccionario con conteos de errores y coordenadas específicas (fila, col)
        SOLO si el error es identificable de forma única (error de un solo bit).
    Args:
        matriz (list[str]): Matriz de bits representando el mensaje.
    Returns:
        dict: Resultado del análisis con conteos y coordenadas.
    """
    # Convierto a int
    mat = [[int(bit) for bit in row] for row in matriz]
    n_rows = len(mat)
    n_cols = len(mat[0])

    error_rows = []
    error_cols = []


    for r in range(n_rows):
        if sum(mat[r]) % 2 != 0:
            error_rows.append(r)

    for c in range(n_cols):
        col_sum = sum(mat[r][c] for r in range(n_rows))
        if col_sum % 2 != 0:
            error_cols.append(c)

    coords = []
    print(f"Error rows {error_rows}")
    print(f"Error cols {error_cols}")
    if len(error_rows) == 1 and len(error_cols) == 1:
        coords.append((error_rows[0], error_cols[0]))


    return {
        "row_error_cont": len(error_rows),
        "col_error_cont": len(error_cols),
        "coords": coords
    }

def correccion(matriz: list[str], resultados: dict):
    """_summary_
        Corrige la matriz SOLO si el error está en el mensaje
        Si el error está en un bit de chequeo, se considera no corregible.
    Args:
        matriz (list[str]): Matriz de bits representando el mensaje.
        resultados (dict): Resultado del análisis con conteos y coordenadas.
    Returns:
        list[str]: Matriz corregida o original si no es corregible.
    """
    error_row_cont = resultados['row_error_cont']
    error_col_cont = resultados['col_error_cont']
    coords = resultados['coords']

    if error_row_cont > 1 or error_col_cont > 1:
        print("  Más de un error detectado. No se corrige.")
        return matriz
    if not coords or error_row_cont == 0 and error_col_cont == 0:
        print("  No se detectaron errores.")
        return matriz

    r, c = coords[0]
    n_cols = len(matriz[0])

    row_paridad = (r == 0)
    col_paridad = (c == n_cols - 1)

    if row_paridad or col_paridad:
        print(f"  Error en un bit de chequeo en ({r}, {c}). No se corrige.")
        return matriz

    caracter = [list(row) for row in matriz]

    bit_actual = caracter[r][c]
    nuevo_bit = '0' if bit_actual == '1' else '1'

    print(f"  Corrigiendo error en posición ({r}, {c}): {bit_actual} -> {nuevo_bit}")
    caracter[r][c] = nuevo_bit

    return ["".join(row) for row in caracter]

ENTRADA = "10111101100111101000001010100101"
SALIDA = "GGHHIIGHGFFGIIFGHGGGGFGGFGGFGFFG"

if __name__ == "__main__":

    matriz_canal, rows, cols = generar_matriz_canal(ENTRADA,SALIDA)
    imprimir_matriz(matriz_canal,rows,cols)
    probs_a_priori = calcular_a_priori(ENTRADA)
    probs_priori = list(probs_a_priori.values())
    print(f"Probabilidades a priori {probs_priori}")

    determinante = es_canal_determinante(matriz_canal)
    print(f"Es determinante: {determinante}")
    sin_ruido = es_canal_sin_ruido(matriz_canal)
    print(f"Es sin ruido: {sin_ruido}")

    equivocacion = calcular_ruido(probs_priori, matriz_canal)
    print(f"Equivocacion: {equivocacion:.4f}")
    perdida = calcular_perdida(probs_priori, matriz_canal)
    print(f"Perdida {perdida:.4f}")
    informacion_mutua = calcular_informacion_mutua(probs_priori, matriz_canal)
    print(f"Info mutua {informacion_mutua:.4f}")
    matriz_reducida = reducir_matriz(matriz_canal)
    filas = [i for i in range(len(matriz_reducida))]
    cols = [j for j in range(len(matriz_reducida[0]))]
    imprimir_matriz(matriz_reducida,filas,cols)

    probs_cap, capacidad = capacidad_max_info(matriz_canal, 0.0001)
    print(f"Probabilidades optimas: {probs_cap}")
    print(f"Capacidad {capacidad:.4f}")
    probabilidad_de_error = prob_error_canal(probs_priori,matriz_canal)
    print(f"Probabilidad de error: {probabilidad_de_error:.4f}")

    matriz_mensaje = [
        ["10111101",
         "10011110",
         "10000010",
         "10100101"]
    ]

    resultados = analizar_paridad(matriz_mensaje)
    matriz_corregida = correccion(matriz_mensaje, resultados)
