
# set hrefs to the page
text_href = '#main-text'
toc_href = '#table-of-contents'
figlist_href = '#figure-list'
trigger_href = '#toc-trigger'
navbar_href = '#navbar'

# responsive web settings
two_column_width = 1200
one_columnn_width = 480
max_youtube_video_width = 500
youtube_video_ratio = 9/16

# declare module variables
text = null
toc = null
figlist = null
trigger = null
navbar = null


# convenience resizing functions with regular API to
# do resizing of elements

set_outer_height = (div, height) ->
  margin = div.outerHeight(true) - div.innerHeight()
  margin += parseInt(div.css('padding-top'))
  margin += parseInt(div.css('padding-bottom'))
  div.height(height - margin)

set_outer_width = (div, width) ->
  margin = div.outerWidth(true) - div.innerWidth()
  margin += parseInt(div.css('padding-left'))
  margin += parseInt(div.css('padding-right'))
  div.width(width - margin)

get_outer_width = (div) -> div.outerWidth(true)

get_spacing_width = (div) -> 
  get_outer_width(div) - get_content_width(div)

get_outer_height = (div) -> div.outerHeight(true)

get_content_width = (div) ->
  width = div.innerWidth()
  width -= parseInt(div.css('padding-left'))
  width -= parseInt(div.css('padding-right'))
  width

get_content_height = (div) ->
  height = div.innerHeight()
  height -= parseInt(div.css('padding-top'))
  height -= parseInt(div.css('padding-bottom'))
  height

get_bottom = (div) -> div.position().top + div.outerHeight(true)

get_right = (div) -> div.position().left + div.outerWidth(true)

get_top = (div) -> div.position().top

get_left = (div) -> div.position().left

set_top = (div, top) -> div.css('top', top)

set_left = (div, left) -> div.css('left', left)

resize_img_dom = (img_dom, width) ->
  img_elem = $(img_dom)
  if img_elem.hasClass('inline-graphic')
    return
  if img_dom.naturalWidth > 0 and img_dom.naturalWidth < width
    img_elem.css('width', '')
  else
    img_elem.css('width', '100%')

get_n_column = () ->
  window_width = $(window).width()
  if window_width <= one_columnn_width
    return 1
  if window_width <= two_column_width
    return 2
  return 3

is_active_trigger_toc = () ->
  trigger.hasClass('active')

toggle_trigger = () ->
  if is_active_trigger_toc()
    trigger.removeClass('active')
  else
    trigger.addClass('active')

get_trigger_toc_width = () ->
  window_width = $(window).width()
  trigger_width = parseInt(1.5*get_outer_width(trigger))
  if get_n_column() == 2
    return Math.round(window_width/2) - trigger_width
  else if get_n_column() == 1
    return window_width - trigger_width

# resize callback for the window
resize_window = () ->
  window_width = $(window).width()
  window_height = $(window).height()
  body = $(document.body)
  body_padding_left = parseInt(body.css('padding-left'))
  body_padding_right = parseInt(body.css('padding-right'))
  navbar_height = get_outer_height(navbar)

  if is_active_trigger_toc()
    toggle_trigger()
  if get_n_column() == 3
    set_left(navbar, body_padding_left)
    set_outer_width(navbar, window_width - body_padding_left - body_padding_right)
  else
    set_left(navbar, 0)
    set_outer_width(navbar, window_width)

  if get_n_column() < 3
    trigger.show()
  else
    trigger.hide()

  height = window_height - navbar_height
  for div in [toc, figlist, text]
    set_top(div, navbar_height)
    set_outer_height(div, height)

  # move toc around
  if get_n_column() < 3
    set_top(toc, 0)
    set_outer_height(toc, window_height)
    width = get_trigger_toc_width()
    set_left(toc, -width)
    set_outer_width(toc, width)
  else if get_n_column() == 3
    toc.width('')
    set_left(toc, body_padding_left)

  # set text and figlist
  if get_n_column() == 1
    figlist.hide()
    set_left(text, 0)
    set_outer_width(text, window_width)
  else if get_n_column() == 2
    figlist.show()
    half_window_width = Math.round(window_width/2)
    set_left(text, 0)
    set_outer_width(text, half_window_width)
    figlist_width = window_width - half_window_width
    set_left(figlist, half_window_width)
    set_outer_width(figlist, figlist_width)
  else # three columns
    figlist.show()
    text.width('')
    left = get_right(toc)
    set_left(text, left)
    left = get_right(text)
    set_left(figlist, left)
    figlist_width = \
        window_width \
        - body_padding_left \
        - body_padding_right \
        - get_outer_width(toc) \
        - get_outer_width(text)
    set_outer_width(figlist, figlist_width)

  # resize images and youtube videos
  if get_n_column() == 1
    $('.fig-in-text iframe[src*="youtube.com"]').each(
      () ->
        iframe = $(this)
        parent_width = iframe.parent().width()
        iframe.width(parent_width)
        height = parent_width * youtube_video_ratio
        iframe.height(height)
    )
  else
    # special youtube video handler
    $('.fig-in-figlist iframe[src*="youtube.com"]').each(() ->
      iframe = $(this)
      parent_width = iframe.parent().width()
      if parent_width < max_youtube_video_width
        iframe.css('width', '100%')
        height = parent_width * youtube_video_ratio
      else
        iframe.css('width', max_youtube_video_width + 'px')
        height = max_youtube_video_width * youtube_video_ratio
      iframe.css('height', height + 'px')
    )

    # and now for images
    for img_dom in $('img')
      img_elem = $(img_dom)
      parent_width = figlist_width - get_spacing_width(img_elem.parent())
      make_resize_fn = (img_dom, parent_width) ->
          () -> resize_img_dom(img_dom, parent_width)
      resize_fn = make_resize_fn(img_dom, parent_width)
      if init == true
        # in case image has not been loaded yet!
        img_elem.load(resize_fn)
      else
        resize_fn(img_dom, parent_width)

trigger_toc = () ->
  toggle_trigger()
  if get_n_column() < 3
    width = get_trigger_toc_width()
    if is_active_trigger_toc()
      move_css = {left: "+="+width}
    else
      move_css = {left: "-="+width}
    for div in [toc, text, navbar, figlist]
      div.animate(move_css, 300)


init = () ->
  text = $(text_href)
  toc = $(toc_href)
  figlist = $(figlist_href)
  trigger = $(trigger_href)
  navbar = $(navbar_href)

  trigger.click(trigger_toc)

  trigger.css('z-index', 1000)
  toc.css('z-index', 700)
  text.css('z-index', 500)
  figlist.css('z-index', 0)

  supplescroll.init_touchscroll()
  supplescroll.build_page(toc_href, text_href, figlist_href)

  $(window).resize(resize_window)
  resize_window()


$(window).ready(init)



