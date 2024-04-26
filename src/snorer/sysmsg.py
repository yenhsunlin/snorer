# Created by Yen-Hsun Lin (Academia Sinica) in 04/2024.
# Copyright (c) 2024 Yen-Hsun Lin.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version (see <http://www.gnu.org/licenses/>).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


class FlagError(Exception):
    
    def __init__(self, message):
        """
        User-defined exception

        Input
        ------
        message: The message to raise when exception ecountered
        """
        self.message = message
        
    def __str__(self):
        return f'{self.message}'
    
    def __call__(self):
        if self.message == 'AFFIDAVIT':
            return f'The truth, the whole truth and nothing but the truth: Boss ALWAYS Dictates!'