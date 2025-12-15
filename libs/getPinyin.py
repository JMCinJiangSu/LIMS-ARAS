#-*- coding:gbk -*-
#from pypinyin import lazy_pinyin
from itertools import chain
from functools import lru_cache
from pypinyin import pinyin, Style

'''
Discription
	
	将中文转化为拼音，便于排序。 
	
'''

@lru_cache(maxsize=128)
def topinyin(instr):
	# 使用pinyin会带上音标，如“中文”返回“zhong1wen2”，改用lazy_pinyin，返回“zhongwen”
	# 2023.06.07-还是改回带音标的，不带的话可能会有“埃克替尼”排在“阿美替尼”前面的情况，而“阿法替尼”会在埃克替尼前，有点奇怪
	return "".join(chain.from_iterable(pinyin(instr, Style.TONE3)))
	#return "".join(chain.from_iterable(lazy_pinyin(instr)))
