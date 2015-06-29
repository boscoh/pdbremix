

# # PageSlab: a module to create
# 1. integrated Table of Contents
# 2. integrated figures
# 3. linked list of figures
# 4. clean resizing of columns
# 5. touchscroll for elems



is_onscreen = (parent_div, div) ->
  x1 = parent_div.offset().top
  x2 = x1 + parent_div.height()
  y1 = div.offset().top 
  y2 = y1 + div.outerHeight(true)
  if x1 <= y1 and y1 <= x2
    return true
  if x1 <= y2 and y2 <= x2
    return true
  if y1 <= x1 and x1 <= y2
    return true
  if y1 <= x2 and x2 <= y2
    return true
  return false


class FigureList
  # text_href must be 'position:relative'
  constructor: (@toc_href, @text_href, @figlist_href) ->
    $.scrollTo.defaults.axis = 'y'
    
    # initialize properties
    @selected_figlink = null
    @selected_header = null
    @selected_headerlink = null
    @is_autodetect_figlink = true
    @is_autodetect_header = true
    @headers = []
    @figlinks = []
    
    $(@text_href).append(
        $('<div>').addClass('page-filler'))

    if @figlist_href != ''
      @transfer_figs()
      @make_reflinks()
      @make_figlinks()
      $(@figlist_href).append(
          $('<div>').addClass('page-filler'))

    if @toc_href != ''
      @make_toc()
      $(@toc_href).append(
          $('<div>').addClass('page-filler'))

    $(@text_href).scroll(() => @scroll_in_text())

    # handle initial hash code
    hash = window.location.hash
    if hash.slice(0, 7) == '#header'
      @select_header($(hash))
    else if hash.slice(0, 4) == '#fig'
      for figlink in @figlinks
        if figlink.attr('href') == hash
          fig = $(hash)
          # wait till all assets have been loaded!
          fig.ready(@select_figlink_fn(figlink))
    else if hash.slice(0, 4) == '#ref'
      for reflink in @reflinks
        if reflink.attr('href') == hash
          ref = $(hash)
          # wait till all assets have been loaded!
          ref.ready(@select_figlink_fn(reflink))

  make_toc: ->
    toc = $(@toc_href)
    div = $('<div>').addClass('toc')
    toc.append(div)

    n_header = 1
    @headers = []
    @headerlinks = {}
    for header_dom in $(@text_href).find('h1, h2, h3, h4')
      # give the header an id to link against
      header = $(header_dom)
      header_id = 'header' + n_header
      n_header += 1
      header.attr('id', header_id)
      @headers.push(header)

      # create a link in the toc
      header_href = '#' + header_id
      headerlink = $('<a>').attr('href', header_href)
      headerlink.append(header.clone().attr('id', ''))
      @headerlinks[header_id] = headerlink
      finish = () =>
        @select_onscreen_figlink_and_figure()
      headerlink.click(@scroll_to_href_in_text_fn(header_href, false, finish))
      div.append(headerlink)

  transfer_figs: () ->
    figlist = $(@figlist_href)
    # in @text_href, move all <div id='fig*'> into @figlist_href
    num_fig = 1
    for div_dom in $(@text_href).find('div')
      div_id = $(div_dom).attr('id')
      if div_id? and div_id[0..2] == 'fig'
        div = $(div_dom)
        # div.prepend('(Figure ' + num_fig + '). ') 
        div.prepend(' ') 
        new_div = div.clone()
        div.addClass('fig-in-text')
        new_div.addClass('fig-in-figlist')
        figlist.append(new_div)
        num_fig += 1

  make_figlinks: () ->
    # in @figlist_href, for <div id='fig*'>, assign label 'Figure.#'
    @i_fig_dict = {}
    @fig_hrefs = []
    @fig_href_from_orig = {}
    @fig_label_dict = {}
    @figlinks = []

    # find all figures in the figlist, and change their id's to figure{n}
    n_fig = 1
    for fig_div_dom in $(@figlist_href).find('div')
      fig = $(fig_div_dom)
      fig_id = fig.attr('id')
      if fig_id? and fig_id[0..2] == 'fig'
        orig_fig_href = '#' + fig_id
        new_fig_href = '#figure' + n_fig
        @fig_href_from_orig[orig_fig_href] = new_fig_href
        @i_fig_dict[new_fig_href] = n_fig
        @fig_hrefs.push(new_fig_href)
        @fig_label_dict[new_fig_href] = $('<span>')
        fig.attr('id', 'figure' + n_fig)
        n_fig += 1

    # find all figlinks, and set their href's and id's
    n_figlink = 1
    for figlink_dom in $(@text_href).find('a[href*="fig"]')
      figlink = $(figlink_dom)

      # make an ID for a figlink so backlinks can point to it
      figlink_id = 'figlink'+n_figlink
      figlink.attr('id', figlink_id)
      figlink.addClass('figlink')
      figlink.click(@select_figlink_fn(figlink))

      orig_fig_href = figlink.attr('href')

      if orig_fig_href of @fig_href_from_orig
        # figure out from the figlink what figure it points
        # to and what the new figure id is
        fig_href = @fig_href_from_orig[orig_fig_href]
        i_fig =  @i_fig_dict[fig_href]
        # figlink_label = 'Figure ' + i_fig + '&rArr;'
        # figlink.html(figlink_label)
        figlink.append('&rArr;')
        figlink.attr('href', fig_href)

        figlink_href = '#'+figlink_id
        reverse_link = $('<a>').append('&lArr;').attr('href', figlink_href)
        select_fig_fn = @select_figlink_fn(figlink)
        finish = ()=>
          select_fig_fn()
          window.location.hash = @selected_figlink.attr('href')
        click_fn = @scroll_to_href_in_text_fn(figlink_href, false, select_fig_fn)
        reverse_link.click(click_fn)

        @figlinks.push(figlink)
        @fig_label_dict[fig_href].append(reverse_link)
        n_figlink += 1
  
    if @figlinks[0]
      @select_figlink(@figlinks[0])
    
    for fig_href, i in @fig_hrefs
      num_fig = i + 1
      fig_label = @fig_label_dict[fig_href]
      $(fig_href).prepend(fig_label)

  make_reflinks: () ->
    @ref_hrefs = []
    @ref_label_dict = {}
    @reflinks = []

    # find all figures in the figlist, and change their id's to figure{n}
    for ref_div_dom in $(@figlist_href).find('a')
      ref = $(ref_div_dom)
      ref_id = ref.attr('id')
      if ref_id? and ref_id[0..3] == 'ref-'
        ref_href = '#' + ref_id
        @ref_hrefs.push(ref_href)
        # initialize DOM object for reverse_links
        @ref_label_dict[ref_href] = $('<span>')

    # find all reflinks, set their href's and id's
    n_reflink = 1
    for reflink_dom in $(@text_href).find('a[href*="ref"]')
      reflink = $(reflink_dom)

      # make an ID for a reflink so backlinks can point to it
      reflink_id = 'reflink'+n_reflink
      reflink.attr('id', reflink_id)
      reflink.append('&rArr;')
      reflink.addClass('reflink')
      reflink.click(@select_figlink_fn(reflink))
      @reflinks.push(reflink)
      n_reflink += 1

      # check for actural ref's pointed to by reflink
      # and makes a dangling DOM object for a reverse_link
      ref_href = reflink.attr('href')
      if ref_href in @ref_hrefs
        reflink_href = '#'+reflink_id
        reverse_link = $('<a>').append('&lArr;').attr('href', reflink_href)
        finish = ()=>
          return
        click_fn = @scroll_to_href_in_text_fn(reflink_href, false, finish)
        reverse_link.click(click_fn)
        @ref_label_dict[ref_href].append(reverse_link)

    for ref_href in @ref_hrefs
      ref = $(@figlist_href).find(ref_href)
      ref_label = @ref_label_dict[ref_href]
      ref.parent().prepend(' ')
      ref.parent().prepend(ref_label)

  select_figlink: (figlink) ->
    if @selected_figlink == figlink
      return
    if @selected_figlink != null
      @selected_figlink.removeClass('active')
      selected_fig_href = @selected_figlink.attr('href')
      $(selected_fig_href).removeClass('active')
    @selected_figlink = figlink
    @selected_figlink.addClass('active')
    selected_fig_href = @selected_figlink.attr('href')
    $(selected_fig_href).addClass('active')

  scroll_to_next_figlink: () ->
    if @is_scrolling_figlist
      # already scrolling so cancel
        return
    finish = () =>
      @is_scrolling_figlist = false
      # we've stopped scrolling now
      if @selected_figlink != @next_figlink
        # but @next_figlink has changed, so do again
        @scroll_to_next_figlink()
    figlist = $(@figlist_href)
    text = $(@text_href)
    if figlist.css('visibility') == 'hidden'
      target = text
    else
      target = figlist
    @is_scrolling_figlist = true
    @select_figlink(@next_figlink)
    fig_href = @selected_figlink.attr('href')
    target.scrollTo(fig_href, 500, finish)

  select_figlink_and_scroll_to_fig: (figlink) ->
    @next_figlink = figlink
    if @selected_figlink == @next_figlink
      return
    @scroll_to_next_figlink()

  select_figlink_fn: (figlink) ->
    (e) => 
      if e? and e.hasOwnProperty('preventDefault')
        e.preventDefault()
      @select_figlink_and_scroll_to_fig(figlink)

  select_header: (header) ->
    @selected_header = header
    header_id = header.attr('id')

    # deselect old header in toc
    if @selected_headerlink != null
      @selected_headerlink.removeClass('active')

    # make header active
    @selected_headerlink = @headerlinks[header_id]
    @selected_headerlink.addClass('active')

    hash = '#' + header_id
    if history.pushState
      history.pushState(null, null, hash)
    else
      window.location.hash = hash

  scroll_to_href_in_text: (href, is_autodetect_figlink, callback) ->
    @is_autodetect_figlink = is_autodetect_figlink
    finish = () => 
      @is_autodetect_figlink=true
      if callback?
        callback()
    delayed_finish = ()->setTimeout(finish, 250)
    settings = { onAfter:delayed_finish, offset:{ top:-15 }}
    $(@text_href).scrollTo(href, 500, settings)

  scroll_to_href_in_text_fn: (href, is_autodetect_figlink, callback) ->
    (e) =>
      e.preventDefault()
      @scroll_to_href_in_text(href, is_autodetect_figlink, callback)

  select_onscreen_figlink_and_figure: () ->
    text = $(@text_href)
    # check if @selected_figlink is onsceen
    if @selected_figlink?
      if is_onscreen(text, @selected_figlink)
        return
    # check if @selected_figlink is onsceen
    onscreen_figlink = null
    for figlink in @figlinks
      if is_onscreen(text, figlink)
        onscreen_figlink = figlink
        break
    if onscreen_figlink?
      @select_figlink_and_scroll_to_fig(onscreen_figlink)

  select_onscreen_header: () ->
    text = $(@text_href)
    # check for onscreen header, and update toc
    # no big changes, so can always run
    onscreen_header = null
    for header in @headers
      if is_onscreen(text, header)
        onscreen_header = header
        break
    if onscreen_header? 
      if @selected_header != onscreen_header
        @select_header(onscreen_header)

  scroll_in_text: () ->
    # $(@text_href) must be position:relative to work
    @select_onscreen_header()
    if @is_autodetect_figlink
      @select_onscreen_figlink_and_figure()



build_page = (toc_href, text_href, figlist_href) ->
    window.figure_list = new FigureList(toc_href, text_href, figlist_href)


# routines to handle the touchscroll on iOS devices

# figure out the DOM element that has triggered
# and scroll it away from the top or bottom of elem
shift_from_edge = (e) ->
  target = e.currentTarget
  bottom = target.scrollTop + target.offsetHeight
  if target.scrollTop == 0
    target.scrollTop = 1
  else if target.scrollHeight == bottom
    target.scrollTop -= 1


init_touchscroll = () ->
  # block whole document from bouncing
  $(document).on(
    'touchmove', 
    (e)->e.preventDefault())
  # allow elements with .touchscroll to bounce
  $('body').on(
      'touchmove', 
      '.touchscroll', 
      (e)->e.stopPropagation())
  # add hack to stop .touchscroll elems from hitting 
  # extremities to avoid triggering the whole page to bounce
  $('body').on(
      'touchstart', 
      '.touchscroll', 
      (e)-> shift_from_edge(e))


# Public API!

window.supplescroll = {
  build_page: build_page,
  init_touchscroll: init_touchscroll
}



