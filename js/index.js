/*
* @Author: XiZhihui
* @Date:   2018-10-31 09:13:49
* @Last Modified by:   XiZhihui
* @Last Modified time: 2018-11-06 17:07:58
* @Description: format content html
*/


function px2rem(px) {
	font_size = getComputedStyle(document.documentElement)["fontSize"]
	font_size = parseFloat(font_size)
	return px / font_size + "rem"
}

/* 
// function Card* format* are used to create html str according template below
<h2 class="Year">year</h2>
<div class="Card">
	<h2 class="Card-title"><a name="">month/category</a></h2>
	<ul class="Card-list">
		<li class="Card-list-item"><a href="url">title</a><span>date</span></li>
		<li class="Card-list-item"><a href="url">title</a><span>date</span></li>
	</ul>
	<div class="Card-more"><button class="Card-more-btn">更多</button></div>
</div>
*/
function CardListItem(card_items, showTime=true) {
	if (!Array.isArray(card_items)) card_items = [card_items]
	card_items.sort((prev, next) => {
		return new Date(next.date).getTime() - new Date(prev.date).getTime()
	})
	if (card_items.length > 10) card_items = card_items.slice(0,10)
	card_items = card_items.map(item => {
		if (!showTime) item.date = ''
		return `<li class="Card-list-item"><a href="${item.url}">${item.title}</a><span>${item.date}</span></li>`
	})
	return card_items.join("")
}

function CardList(card_items) {
	card_items = CardListItem(card_items)
	return `<ul class="Card-list js-card-list">${card_items}</ul>`
}

function CardTitle(title, year) {
	let name = title
	if (year) {
		name = year + title
		title = "January February March April May June July August September October November December".split(" ")[parseInt(title-1)]
	}
	return `<h2 class="Card-title"><a name="${name}">${title}</a></h2>`
}

function Card(data, type, year) {
	let titles = Object.keys(data)
	titles.sort().reverse()

	let card_str = ''
	let cards = []
	let moreOrLess = ''

	titles.forEach(title => {
		let items = data[title]
		let card_list = CardList(items)
		title_str = CardTitle(title, year)

		if (type != 'about' && items.length > 10) {
			moreOrLess = `<div class="Card-more" data-index="0" data-year="${year}" data-title="${title}" data-type="${type}">
						<button class="js-more-btn" onclick="CardMoreOrLess(event)">更多</button>
						<button class="js-less-btn" style="display: none" onclick="CardMoreOrLess(event)">收起</button>
					</div>`
		}
		card_str = `<div class="Card">
					${title_str}
					${card_list}
					${moreOrLess}
				</div>`
		cards.push(card_str)
	})
	// cards.reverse()
	return cards.join("")
}
	// 给按钮定义响应函数
function CardMoreOrLess(evt) {
	if (!evt) evt = window.event
	// 如果不是按钮直接返回
	let target = evt.target || evt.srcElement

	let txt = target.innerText
	let parent = target.parentNode
	let data = parent.dataset
	let card_list = parent.previousElementSibling
	let height
	if ("height" in parent.dataset) {
		height = parent.dataset.height
	} else {
		height = parseFloat(getComputedStyle(card_list,"")["height"])
		parent.dataset.height = height
	}

	index = parseInt(data.index)

	// 得到当前card的数据, current_json 最后是个数组
	// current_json = window[data.type + "_json"]
	// if (data.year == "undefined") {
	// 	current_json = current_json[data.title]
	// } else {
	// 	current_json = current_json[data.year][data.title]
	// }
	current_json = get_json(total_json)[data.type]
	if (data.type == "index") {
		current_json = current_json[data.year][data.title]
	}

	if (txt == "更多") {
		// 根据 data.index 确定已经显示了多少,再显示10个
		let items = current_json.slice((index+1)*10, (index+2)*10)
		let targetHeight = height / 10 * ((index + 1)* 10 + items.length)

		card_list.innerHTML += CardListItem(items)
		parent.dataset.index = index+1

		target.nextElementSibling.style.display = 'inline-block'

		// 添加过渡效果
		card_list.style.height = targetHeight + "px"

		// 需要的数据没有了，不显示更多按钮
		if ((index+2)*10 > current_json.length) {
			target.style.display = "none"
		}
	} else {
		parent.dataset.index = 0
		target.previousElementSibling.style.display = 'inline-block'
		target.style.display = "none"
		card_list.style.height = height + "px"
		setTimeout(()=> {
			card_list.innerHTML = CardListItem(current_json)
		}, 1000)
	}
	target.focus()
}
function format_content(data, type) {
	let html_str = ''
	if (type == "index") {
		years = Object.keys(data)
		years.sort().reverse()
		years.forEach(year => {
			html_str += `<h2 class="Year"><a name="${year}">${year}</a></h2>`
			html_str += Card(data[year], type, year)
		})
	} else {
		html_str += Card(data, type)
	}
	return html_str
}
//——————————————————————————————————————————————————
// 目录html
function format_table(data, type) {
	let html_str = []
	let months = "January February March April May June July August September October November December".split(" ")
	if (type == "index") {
		let temp = [],
			temp_str = ''
		for (year in data) {
			temp_str = `<h3 class="Table-title"><a href="#${year}">${year}</h3>`
			for (month in data[year]) {
				month_str = `<h4 class="Table-title-item"><a href="#${year+month}">${months[parseInt(month)-1]} (${data[year][month].length})</a></h4>`
				temp.push(month_str)
			}
			temp.reverse() // 最近的内容出现页面上方
			temp_str += temp.join("")
			html_str.push(temp_str)
		}
	}else {
		for (cate in data) {
			let cate_num = cate
			if (type != "about") cate_num = `${cate} (${data[cate].length})`
			cate_str = `<h4 class="Table-title-item"><a href="#${cate}">${cate_num}</a></h4>`
			html_str.push(cate_str)
		}
	}
	return html_str.join("")
}

function format_article_table(tables) {
	// 一个小结的结构：
	// h2(h2)
	// 	 h3(ul )
	//     h4 (li)
	//   如果出现 h3,就另开一个
	// 如果出现 h2 就另开一个这样的结构
	// [['h2', 'gokegg', 'GO/KEGG富集分析'], 
		// ['h3', '1ensembl-identrezid', '1.对Ensembl ID进行转换,得到对应的基因名和EntrezID'], 
		// ['h3', '2', '2. 进行注释'], 
		// ['h4', '23-gsea', '2.3 GSEA富集分析'], 
		// ['h3', '3', '3.进行可视化']]
	current = []
	totals = []
	html = ''
	for (let idx = 0, tablen = tables.length; idx < tablen; ) {
		table = tables[idx]
		idx++
		if (table[0] == 'h2') {
			totals.push(current)
			current = [table]
		} else {
			current.push(table)
		}
	}
	totals.push(current)
	console.log(totals)
	while (totals.length > 0) {
		current = totals.shift()
		temp_h3 = ''
		temp_h4 = ''

		while (current.length > 0) {
			let temp = current.pop()
			if (temp[0] == 'h4') {
				temp_h4 = `<li><a href='#${temp[1]}'>${temp[2]}</a></li>` + temp_h4
			} else if (temp[0] == 'h3') {
				if (temp_h4.length) temp_h4 = `<ul>${temp_h4}</ul>`
				temp_h3 = `<li><h4><a href='#${temp[1]}'>${temp[2]}</a></h4>${temp_h4}</li>` + temp_h3
				temp_h4 = ''
			} else {
				if (temp_h3.length) temp_h3 = `<ul>${temp_h3}</ul>`
				html += `<div><h3><a href='#${temp[1]}'>${temp[2]}</a></h3>${temp_h3}</div>`
			}
		}
	}

	return html
}
//———————————————————————————————————————————————————
// 页面切换
function get_json(raw_json) {
	index_json = {}
	category_json = {}
	raw_json.forEach(item => {
		// index_json
		[year, month, day] = item.date.split("/")
		if (year in index_json) {
			if (month in index_json[year]) {
				index_json[year][month].push(item)
			} else {
				index_json[year][month] = [item]
			}
		} else {
			index_json[year] = {}
		}
	})
	return {
		index: index_json
	}
}

function change_page(evt) {
	let target = evt.target || evt.srcElement
	if (target.tagName != "LI") return

	let pages = Array.from(target.parentNode.querySelectorAll("li"))
	// 页面内容和左侧目录
	let content = '', nav = '', current_json

	// 消除当前页的 js 钩子
	pages.forEach(page => {
		page.classList.remove("js-current-page")
	})
	// 为新当前页添加 js 钩子
	target.classList.add("js-current-page")
	txt = target.innerText
	location.hash = txt

	document.querySelector(".Search").style.display = "block"

	jsons = get_json(total_json)
	// 根据新的当前页更新页面
	switch (txt) {
		case "Timeline":
			timeline_json = jsons.index
			content = format_content(timeline_json, "index")
			nav = format_table(timeline_json, "index")
			break
		case "Category":
			content = format_content(category_json, "category")
			nav = format_table(category_json, "category")
			break
		case "Tags":
			content = format_content(tags_json, "tags")
			nav = format_table(tags_json, "tags")
			break
		case "About":
			content = format_content(about_json, "about")
			nav = format_table(about_json, "about")
			// 清除搜索框
			document.querySelector(".Search").style.display = "none"
			break
	}
	document.querySelector("content#content").innerHTML = content
	document.querySelector("nav#table").innerHTML = nav
	// 为点击更多的动画做准备
	Array.from(document.querySelectorAll(".js-card-list")).forEach(card => {

		height = parseFloat(getComputedStyle(card)["height"])
		card.style.height = px2rem(height)
	})
}
//———————————————————————————————————————————————————
//响应搜索
function search(evt) {
	if (!evt) evt = window.event

	let target = evt.target || evt.srcElement
	let aim = target.value.toLowerCase()
	let current = document.querySelector(".js-current-page").innerText
	let filtered = []
	let result_list = document.querySelector("main div.Search ul.Search-result")

	if (current == "Timeline") {
		for (year in index_json) {
			let data = index_json[year]
			for (key in data) {
				let temp = data[key].filter( d => {
					return d.title.toLowerCase().indexOf(aim) > -1 || d.date.toLowerCase().indexOf(aim) > -1
				})
				filtered = filtered.concat(temp)
			}
		}
	} else {
		let data = window[current.toLowerCase() + "_json"]
		for (key in data) {
			let temp = data[key].filter( d => {
				return d.title.toLowerCase().indexOf(aim) > -1 || d.date.toLowerCase().indexOf(aim) > -1
			})
			filtered = filtered.concat(temp)
		}
	}

	if (!aim) {
		result_list.innerHTML = ""
		result_list.classList.remove("show-result")
	} else if (filtered.length == 0) {
		result_list.classList.add("show-result")
		result_list.innerHTML = '<li class="Card-list-item">没有找到相关文章</li>'
	} else {
		result_list.innerHTML = CardListItem(filtered, false)
		result_list.classList.add("show-result")
	}
}

//———————————————————————————————————————————————————
//调整某些文章页面的导航
function resize_nav(ele) {
	let maxheight = parseFloat(window.innerHeight)
	let height = parseFloat(getComputedStyle(ele)["height"])
	if (height > maxheight * 0.6) {
		ele.style.zoom = '0.4'
	}
}

// 响应式：根据 font-size 和 rem 调节px
function adaption(doc, win, designWidth) {
	designWidth = designWidth || 1903
    var docE1 = doc.documentElement,
        resizeEvt = 'orientationchange' in window ? 'orientationchange' : 'resize',
        recalc = function(){
            var clientWidth = docE1.clientWidth;
            if(!clientWidth) return;
            docE1.style.fontSize = 20 * (clientWidth / designWidth) + 'px';
            console.log(20 * (clientWidth / designWidth) + 'px')       
        };

    if (!doc.addEventListener) return;
    win.addEventListener(resizeEvt,recalc,false);
    doc.addEventListener('DOMContentLoaded',(evt) => {
    	recalc()
    	// 左侧导航resize
    	resize_nav(document.querySelector("main #table"))
    }, false);
}

// 初始化
function page_init() {
	// 监听滚动事件
	document.addEventListener("scroll", (evt) => {
		let nav_table = document.getElementById("table")
		if (!nav_table) {
			document.removeEventListener("scroll")
		}
		let motto = document.querySelector(".Motto"),
			scroll_threshold = getComputedStyle(motto)["height"]
		let header = document.querySelector(".Header"),
			top = getComputedStyle(header)["height"]
		let position = getComputedStyle(nav_table)["position"]

		if (parseFloat(document.documentElement.scrollTop) >= parseFloat(scroll_threshold)) {
			if (nav_table.style.position == "fixed") return
			nav_table.style.position = "fixed"
			nav_table.style.top = parseFloat(top) + 10 + "px"
		} else {
			nav_table.style = "position: abosulte"
		}
	}, false)

	// // 响应式
	// adaption(document, window)
}

page_init()