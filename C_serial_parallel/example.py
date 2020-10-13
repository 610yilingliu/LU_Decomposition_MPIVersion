to_find = "WARC-Target-URI:"
shift = len(to_find)
html_dict = dict()
# block 是那一大坨内容的字符串
for block in response:
    start = block.find(to_find)
    # http 链接起始位置.前面那个start是W的位置
    start += shift
    # 初始化一个空字符串
    url = ""
    # 你直接用start也可以,但是cur(current)比较方便理解.表示当前在读哪个位置
    cur = start
    # 如果还没读到换行,就是url没结束. and是以防万一比如有什么奇葩情况(targeturl是那个block的最后一行)
    while block[cur] != '\n' and cur < len(block):
        # url 中加上当前读到的字符
        url = url + block[cur]
        # 指针前进一位
        cur += 1
    # 读内容...假设你后面读的内容叫html
    html_dict[url] = html
