import numpy as np

DEFAULT_LABEL_FONT_SIZE: str = 7

def to_counter_clockwise(clockwise_angle):
    counterclockwise_angle = (90 - clockwise_angle) % 360
    return counterclockwise_angle

def circular_midpoint(start, end, total_count):
    clockwise_distance = (end - start) % total_count
    counterclockwise_distance = (start - end) % total_count
    
    if clockwise_distance <= counterclockwise_distance:
        midpoint = (start + clockwise_distance / 2) % total_count
    else:
        midpoint = (start - counterclockwise_distance / 2) % total_count
        
    return midpoint

def circular_length(start, end, total_count):
    clockwise_distance = (end - start) % total_count
    counterclockwise_distance = (start - end) % total_count
    
    length = min(clockwise_distance, counterclockwise_distance)
    
    return length

