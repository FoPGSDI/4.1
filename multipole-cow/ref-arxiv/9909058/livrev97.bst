% BibTeX Living Reviews bibliography style `livrev'
%
% Version 1.1 for BibTeX versions 0.99a or later, LaTeX version 2.09.
%
% Lee Ashton Wild  10/9/97.
%
% Supported formats : article.
%

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
    onlinetitle
    onlineaddress
    onlinemonth
    onlineyear
    cited
    fileformat
    keywords
  }
  {}
  { label }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
	{ add.period$ write$
	  newline$
	  "\newblock " write$
	}
	{ output.state before.all =
	    'write$
	    { add.period$ " " * write$ }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" write$
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
	'skip$
	{ after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

FUNCTION {bold}
{ duplicate$ empty$
    { pop$ "" }
    { "{\bf " swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{vv~}{ll}{, jj}{, ff}" format.name$ 't :=
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
	    { numnames #1 >
		{ "," * }
		'skip$
	      if$
	      t "others" =
		{ " et~al." * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.onlinenames}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{ff~}{vv~}{ll}{, jj}" format.name$ 't :=
      nameptr #1 >
	{ namesleft #1 =
	    { ", et al." * }
	    { "" * }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { "\comment{author}" author format.names * }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { "\comment{editor}" editor format.names *
      editor num.names$ #1 >
	{ ", eds." * }
	{ ", ed." * }
      if$
    }
  if$
}

FUNCTION {format.onlineeditors}
{ editor empty$
    { "" }
    { "For a related online version see: \comment{editor}" editor format.onlinenames *
      editor num.names$ #1 >
	{ ", eds." * }
	{ ", ed." * }
      if$
    }
  if$
}

FUNCTION {format.onlineauthors}
{ author empty$
    { editor empty$
        { "For a related online version see: "
           "No author or editor in " cite$ * warning$ }
        { format.onlineeditors }
      if$
    }
    { "For a related online version see: \comment{author}" author format.onlinenames * }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { "\comment{title}``" title * "''" *}
  if$
}

FUNCTION {format.onlinetitle}
{ onlinetitle empty$
    { "" }
    { "\comment{onlinetitle}``" onlinetitle * "''" *}
  if$
}

FUNCTION {format.onlineaddress}
{ onlineaddress empty$
    { "" }
    { "\comment{onlineaddress}" onlineaddress *}
  if$
}

FUNCTION {format.address}
{ address empty$
    { "" }
    { "\comment{address}" address *}
  if$
}

FUNCTION {format.journal}
{ journal empty$
    { "" }
    { "\comment{journal}" journal emphasize * }
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {format.date}
{ year empty$
    { month empty$
	{ "" }
	{ "there's a month but no year in " cite$ * warning$
	  month
	}
      if$
    }
    { month empty$
	{ "(\comment{year}" year * ")" * }
	{ "(\comment{month}" month * ", \comment{year}" * year * ")" * }
      if$
    }
  if$
}

FUNCTION {format.onlinedate}
{ onlineyear empty$
    { onlinemonth empty$
	{ "" }
	{ "there's a month but no year in " cite$ * warning$
	  onlinemonth
	}
      if$
    }
    { onlinemonth empty$
	{ "(\comment{onlineyear}" onlineyear * ")" * }
	{ "(\comment{onlinemonth}" onlinemonth * ", \comment{onlineyear}" * onlineyear * ")" * }
      if$
    }
  if$
}

FUNCTION {format.btitle}
{ "\comment{title}{\em " title * "}" *
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.keywords}
{ keywords empty$
   { new.sentence "\keywords{}" }
   { new.sentence "\keywords{" keywords * "}" * }
   if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "volume" volume tie.or.space.connect
      series empty$
	'skip$
	{ " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
	{ series field.or.null }
	{ output.state mid.sentence =
	    { "number" }
	    { "Number" }
	  if$
	  number tie.or.space.connect
	  series empty$
	    { "there's a number but no series in " cite$ * warning$ }
	    { " in " * series * }
	  if$
	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
	{ edition "l" change.case$ " edition" * }
	{ edition "t" change.case$ " edition" * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
	{ "\comment{pages}" pages n.dashify tie.or.space.connect }
	{ "\comment{pages}" pages tie.or.space.connect }
      if$
    }
  if$
}

FUNCTION {format.fileformat.cited}
{ fileformat empty$
     { "there's no fileformat in " cite$ * warning$ }
     { cited empty$
         { "there's no cited in " cite$ * warning$ }
         { "\comment{fileformat}[" fileformat * "]: cited on \comment{cited}" * cited * }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume empty$
    { "there is no volume in " cite$ * warning$ }
    { "\comment{volume}" volume bold * }
  if$
  number empty$
    'skip$
    { "\comment{number}(" number * ")" * *
      volume empty$
	{ "there's a number but no volume in " cite$ * warning$ }
	'skip$
      if$
    }
  if$
  pages empty$
    { "there are no pages in " cite$ * warning$ }
    { duplicate$ empty$
        { pop$ format.pages }
         { ", \comment{pages}" * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.publisher.address.year}
{ publisher empty$
     { address empty$
         { year empty$
            { "there is no publisher, address and year in " cite$ * warning$}
            { "there is no publisher and address in " cite$ * warning$
               "(\comment{year}" year * ")" *}
            if$
         }
         { year empty$
             { "(\comment{address}" address * ")" *
                "there's an address but no publisher and year in " cite$ * warning$ }
             { "there's no publisher in " cite$ * warning$
               "(\comment{address}" address * ", \comment{year}" * year * ")" *
             }
           if$
         }
         if$
      }
     { address empty$
         { year empty$
             { "there's no address and year in " cite$ * warning$ 
                "(\comment{publisher}" publisher * ")" *
             }
             { "there's no address in " cite$ * warning$
                "(\comment{publisher}" publisher * ", \comment{year}" * year * ")" *
             }
            if$
          }
          { year empty$
               { "there's no year in " cite$ * warning$ 
                  "(\comment{publisher}" publisher * ", \comment{address}" * address * ")" * }
               { "(\comment{publisher}" publisher * ", \comment{address}" * address * ", \comment{year}" * year * ")" * }
            if$
         }
      if$
    }
  if$
}

FUNCTION {format.school.address.year}
{ school empty$
     { address empty$
         { year empty$
            { "there is no school, address and year in " cite$ * warning$}
            { "there is no school and address in " cite$ * warning$
               "(\comment{year}" year * ")" *}
            if$
         }
         { year empty$
             { "(\comment{address}" address * ")" *
                "there's an address but no school and year in " cite$ * warning$ }
             { "there's no school in " cite$ * warning$
               "(\comment{address}" address * ", \comment{year}" * year * ")" *
             }
           if$
         }
         if$
      }
     { address empty$
         { year empty$
             { "there's no address and year in " cite$ * warning$ 
                "(\comment{school}" school * ")" *
             }
             { "there's no address in " cite$ * warning$
                "(\comment{school}" school * ", \comment{year}" * year * ")" *
             }
            if$
          }
          { year empty$
               { "there's no year in " cite$ * warning$ 
                  "(\comment{school}" school * ", \comment{address}" * address * ")" * }
               { "(\comment{school}" school * ", \comment{address}" * address * ", \comment{year}" * year * ")" * }
            if$
         }
      if$
    }
  if$
}

FUNCTION {format.tr.institution.address.year}
{ institution empty$
     { address empty$
         { year empty$
            { "there is no institution, address and year in " cite$ * warning$}
            { "there is no institution and address in " cite$ * warning$
               "(\comment{year}" year * ")" *}
            if$
         }
         { year empty$
             { "(\comment{address}" address * ")" *
                "there's an address but no institution and year in " cite$ * warning$ }
             { "there's no institution in " cite$ * warning$
               "(\comment{address}" address * ", \comment{year}" * year * ")" *
             }
           if$
         }
         if$
      }
     { address empty$
         { year empty$
             { "there's no address and year in " cite$ * warning$ 
                "(\comment{institution}" institution * ")" *
             }
             { "there's no address in " cite$ * warning$
                "(\comment{institution}" institution * ", \comment{year}" * year * ")" *
             }
            if$
          }
          { year empty$
               { "there's no year in " cite$ * warning$ 
                  "(\comment{institution}" institution * ", \comment{address}" * address * ")" * }
               { "(\comment{institution}" institution * ", \comment{address}" * address * ", \comment{year}" * year * ")" * }
            if$
         }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
	{ "chapter" }
	{ type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
	'skip$
	{ ", \comment{pages}" * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
	{ "in \comment{booktitle}" booktitle emphasize * }
	{ "in " format.editors * ", \comment{booktitle}" * booktitle emphasize * }
      if$
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    'skip$
    'type
  if$
  number empty$
    'skip$
    { "\comment{number}Report no. " number * }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
	{ "need key or journal for " cite$ * " to crossref " * crossref *
	  warning$
	  ""
	}
	{ "In \comment{journal}{\em " journal * "\/}" * }
      if$
    }
    { "in \cite{???}" key * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
	'skip$
	{ editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " et~al." * }
	    { " and " * editor #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Volume" volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
	    { "\comment{series}{\em " * series * "\/}" * }
	  if$
	}
	{ key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ booktitle empty$
	    { "need editor, key, or booktitle for " cite$ * " to crossref " *
	      crossref * warning$
	      ""
	    }
	    { "In \comment{booktitle}{\em " booktitle * "\/}" * }
	  if$
	}
	{ "in \comment{??}" key * }
      if$
    }
    { "in \comment{editor}" format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {article}
{ output.bibitem
   "\comment{article}" write$ newline$ ""
   format.authors "author" output.check
   format.title "title" output.check
   crossref missing$
    { format.journal "journal" output.check
       format.vol.num.pages "pages" output.check
       format.date "year" output.check
       new.sentence
       onlineaddress empty$
           'skip$
           { format.onlineauthors "author" output.check
              format.onlinetitle "onlinetitle" output.check
              format.onlinedate "year" output.check
              format.fileformat.cited "cited" output.check
              format.onlineaddress "onlineaddress" output.check
              new.sentence
            }
         if$
    }
    { format.article.crossref output.nonnull
       format.pages output
       new.sentence 
    }
  if$
  note output
  format.keywords output
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  "\comment{book}" write$ newline$ ""
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      format.number.series output
      format.publisher.address.year "year" output.check
    }
    { new.sentence
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  new.sentence
  onlineaddress empty$
      'skip$
      { format.onlineauthors "author" output.check
         format.onlinetitle "onlinetitle" output.check
         format.onlinedate "year" output.check
         format.fileformat.cited "cited" output.check
         format.onlineaddress "onlineaddress" output.check
       }
    if$
  note output
  format.keywords output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  "\comment{inbook}" write$ newline$ ""
  format.authors "author" output.check
  format.title "title" output.check
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
       format.bvolume output
       format.chapter.pages "chapter and pages" output.check
       format.number.series output
       format.publisher.address.year "year" output.check
       format.edition output
       new.sentence
       onlineaddress empty$
           'skip$
           { format.onlineauthors "author" output.check
              format.onlinetitle "onlinetitle" output.check
              format.onlinedate "year" output.check
              format.fileformat.cited "cited" output.check
              format.onlineaddress "onlineaddress" output.check
            }
         if$
    }
    { format.chapter.pages "chapter and pages" output.check
       format.book.crossref output.nonnull
       format.edition output 
    }
  if$
  new.sentence
  note output
  format.keywords output 
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  "\comment{inproceedings}" write$
  newline$
  ""
  format.authors "author" output.check
  format.title "title" output.check
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      format.pages output
      address empty$
        { organization publisher new.sentence.checkb
           organization output
           publisher output
           format.date "year" output.check
        }
        {
           format.publisher.address.year output
           organization output
        }
      if$
      new.sentence
      onlineaddress empty$
          'skip$
          { format.onlineauthors "author" output.check
             format.onlinetitle "onlinetitle" output.check
             format.onlinedate "year" output.check
             format.fileformat.cited "cited" output.check
             format.onlineaddress "onlineaddress" output.check
           }
        if$
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
  new.sentence
  note output
  format.keywords output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {lr_url}
{ output.bibitem
  "\comment{lr_url}" write$
  newline$
  ""
  author empty$
    { "\comment{organization}" write$ organization "organization" output.check }
    { "\comment{author}" write$ author "author" output.check }
   if$
  format.onlinetitle "onlinetitle" output.check
  format.onlinedate "onlineyear" output.check 
  format.fileformat.cited "cited" output.check
  format.onlineaddress "onlineaddress" output.check
  new.sentence
  note output
  format.keywords output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  "\comment{mastersthesis}" write$
  newline$
  ""
  format.authors "author" output.check
  format.btitle "title" output.check
  "Master's thesis" format.thesis.type output.nonnull
  format.school.address.year "year" output.check
  new.sentence
  note output
  onlineaddress empty$
      'skip$
      { new.sentence
         format.onlineauthors "author" output.check
         format.onlinetitle "onlinetitle" output.check
         format.onlinedate "year" output.check
         format.fileformat.cited "cited" output.check
         format.onlineaddress "onlineaddress" output.check
       }
    if$
  format.keywords output 
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  "\comment{misc}" write$ newline$ ""
  format.authors output
  format.title output
  howpublished output
  note output
  format.date output
  format.keywords output
  fin.entry
  empty.misc.check
}

FUNCTION {online}
{ output.bibitem
  "\comment{online}" write$ newline$ ""
  format.authors "author" output.check
  format.onlinetitle "onlinetitle" output.check
  format.onlinedate "onlineyear" output.check
  format.fileformat.cited "cited" output.check
  format.onlineaddress "onlineaddress" output.check
  new.sentence
  note output
  format.keywords output
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  "\comment{phdthesis}" write$ newline$ ""
  format.authors "author" output.check
  format.btitle "title" output.check
  "PhD thesis" format.thesis.type output.nonnull
  format.school.address.year "year" output.check
  new.sentence
  note output
  onlineaddress empty$
      'skip$
      { new.sentence
         format.onlineauthors "author" output.check
         format.onlinetitle "onlinetitle" output.check
         format.onlinedate "year" output.check
         format.fileformat.cited "cited" output.check
         format.onlineaddress "onlineaddress" output.check
       }
    if$
  format.keywords output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  "\comment{proceedings}" write$
  newline$
  ""
  editor empty$
    { organization output }
    { format.editors output.nonnull }
  if$
  format.btitle "title" output.check
  format.bvolume output
  format.number.series output
  address empty$
    { editor empty$
       { publisher new.sentence.checka }
       { organization publisher new.sentence.checkb
          organization output
       }
      if$
      publisher output
      format.date "year" output.check
    }
    { format.publisher.address.year output
      editor empty$
         'skip$
         { organization output }
      if$
    }
  if$
  new.sentence
  onlineaddress empty$
      'skip$
      { format.onlineauthors "author" output.check
         format.onlinetitle "onlinetitle" output.check
         format.onlinedate "year" output.check
         format.fileformat.cited "cited" output.check
         format.onlineaddress "onlineaddress" output.check
         new.sentence 
       }
    if$ 
  note output
  format.keywords output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  "\comment{techreport}" write$ newline$ ""
  format.authors output
  format.btitle "title" output.check 
  format.tr.institution.address.year output
  format.tr.number output.nonnull 
  new.sentence
  onlineaddress empty$
      'skip$
      { format.onlineauthors "author" output.check
         format.onlinetitle "onlinetitle" output.check
         format.onlinedate "year" output.check
         format.fileformat.cited "cited" output.check
         format.onlineaddress "onlineaddress" output.check
         new.sentence 
       }
    if$ 
  note output
  format.keywords output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  "\comment{unpublished}" write$
  newline$
  ""
  format.authors "author" output.check
  format.title "title" output.check
  note "note" output.check
  format.date output
  format.keywords output
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}

MACRO {feb} {"February"}

MACRO {mar} {"March"}

MACRO {apr} {"April"}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"August"}

MACRO {sep} {"September"}

MACRO {oct} {"October"}

MACRO {nov} {"November"}

MACRO {dec} {"December"}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { nameptr #1 >
	{ "   " * }
	'skip$
      if$
      s nameptr "{vv{ } }{ll{ }}{  ff{ }}{  jj{ }}" format.name$ 't :=
      nameptr numnames = t "others" = and
	{ "et al" * }
	{ t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
	{ "to sort, need author or key in " cite$ * warning$
	  ""
	}
	{ key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
	{ key empty$
	    { "to sort, need author, editor, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.organization.sort}
{ author empty$
    { organization empty$
	{ key empty$
	    { "to sort, need author, organization, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ "The " #4 organization chop.word sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.organization.sort}
{ editor empty$
    { organization empty$
	{ key empty$
	    { "to sort, need editor, organization, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ "The " #4 organization chop.word sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {presort}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
	'editor.organization.sort
	{ type$ "manual" =
	    'author.organization.sort
	    'author.sort
	  if$
	}
      if$
    }
  if$
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT

STRINGS { longest.label }

INTEGERS { number.label longest.label.width }

FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}

FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}

EXECUTE {initialize.longest.label}

ITERATE {longest.label.pass}

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{"  longest.label  * "}" * write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}

