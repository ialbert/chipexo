import re

ROMAN = ['0', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X',
         'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII', 'XVIII', 'XIX', 'XX',
         'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI', 'XXVII', 'XXVIII', 'XXIX', 'XXX']



def zeropad_to_numeric(data):
    return re.sub(r'chr0(\d)', r'chr\1', data)
    
def numeric_to_zeropad(data):
    return re.sub(r'chr(\d[^\d])', r'chr0\1', data)
    
def roman_to_numeric(data):
    for num, roman in enumerate(ROMAN):
        data = re.sub('chr%s([^IVX])' % roman, r'chr%d\1' % num, data)
    return data

def numeric_to_roman(data):
    for num, roman in enumerate(ROMAN):
        data = re.sub('chr%d([^\d])' % num, r'chr%s\1' % roman,  data)
    return data
    
    
if __name__ == '__main__':
    print zeropad_to_numeric('asdkasjldaj chr02 asdads chr14 asdsahdskj chr3 zxczxc')
    print numeric_to_zeropad('asdkasjldaj chr02 asdads chr14 asdsahdskj chr3 zxczxc')
    print numeric_to_roman('asfdhfkj chr2 asdjhaskj chr14 adshakdjh chr8 asdas')
    print roman_to_numeric('asdasd chrIV jdhakdjsah chrXI adasjhd chrII asd')