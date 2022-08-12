def get_hms(timeinsec: int) -> list:
    '''
    Split time in seconds into hours, minutes and seconds
    '''

    hours = int( timeinsec / 3600 )
    minutes = int( (timeinsec - hours * 3600) / 60 )
    seconds = int( timeinsec - hours * 3600 - minutes * 60 )

    return([hours, minutes, seconds])

hours, minutes, seconds = get_hms(22611)

print(hours, minutes, seconds)