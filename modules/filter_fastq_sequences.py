"""
Модуль для работы с FASTQ последовательностями
"""


def calculate_gc_content(sequence: str) -> float:
    """
    Вычисляет GC состав последовательности в процентах
    
    Args:
        sequence: Нуклеотидная последовательность
    
    Returns:
        GC состав в процентах
    """
    if not sequence:
        return 0.0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100


def calculate_quality_score(quality_string: str) -> float:
    """
    Вычисляет среднее качество последовательности в phred33
    
    Args:
        quality_string: Строка качества в формате phred33
    
    Returns:
        Среднее качество последовательности
    """
    if not quality_string:
        return 0.0
    
    total_quality = sum(ord(char) - 33 for char in quality_string)
    return total_quality / len(quality_string)


def parse_bounds(bounds, default_low: float = 0):
    """
    Парсит границы фильтрации
    
    Args:
        bounds: Интервал или верхняя граница
        default_low: Значение по умолчанию для нижней границы
    
    Returns:
        Кортеж (нижняя_граница, верхняя_граница)
    """
    if isinstance(bounds, (int, float)):
        return (default_low, float(bounds))
    elif isinstance(bounds, tuple) and len(bounds) == 2:
        return (float(bounds[0]), float(bounds[1]))
    else:
        raise ValueError("Границы должны быть числом или кортежем из двух чисел")


def filter_fastq_sequences(
    seqs: dict,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold: int = 0
) -> dict:
    """
    Фильтрует FASTQ последовательности по заданным критериям
    
    Args:
        seqs: Словарь FASTQ последовательностей
        gc_bounds: Границы GC состава
        length_bounds: Границы длины
        quality_threshold: Порог качества
    
    Returns:
        Отфильтрованный словарь последовательностей
    """
    gc_min, gc_max = parse_bounds(gc_bounds)
    length_min, length_max = parse_bounds(length_bounds, default_low=0)
    
    filtered_seqs = {}
    
    for seq_name, (sequence, quality) in seqs.items():
        # Проверка длины
        seq_length = len(sequence)
        if not (length_min <= seq_length <= length_max):
            continue
        
        # Проверка GC состава
        gc_content = calculate_gc_content(sequence)
        if not (gc_min <= gc_content <= gc_max):
            continue
        
        # Проверка качества
        avg_quality = calculate_quality_score(quality)
        if avg_quality < quality_threshold:
            continue
        
        filtered_seqs[seq_name] = (sequence, quality)
    
    return filtered_seqs
