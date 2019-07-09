%% File: ascelike.bst
%%
%% Version 2.2, February 14, 2013
%%
%%
%% An unofficial bibliography style-file for ASCE-like bibliographies.
%%   Adapted from chicago.bst (24-jan-88 version) of Glenn Paulley
%%   This bst style should be used with document class ascelike.cls 
%%   for citations that come out like (Jones, 1986) or variations of such.
%%
%% For use with the BibTeX bibliographic formatting tool, Copyright (C) 1988, 
%%   all rights reserved.  Bibtex is a program for formatting bibliographic
%%   entries.  It uses an author's *.bib file which contains a data base of
%%   information on each entry.  The entries are formatted in a style
%%   specified in a *.bst file, such as this file, so that the editorial
%%   style of the particular book or journal publisher is met.  For example,
%%   this ascelike.bst file prepares bibliographies in an ASCE-like style, 
%%   with quotation marks around titles, etc.
%%
%% To use this ascelike.bst file, you must prepare a *.bib file of your
%%   bibliographic data base.  See The Latex book "LaTeX:  A Document 
%%   Preparation System" by Leslie Lamport, or in Oren Patashnik's "BibTEXing"
%%   document (contained in the btxdoc.* files on CTAN and in online web
%%   versions, search for "btxdox" or "bibtexing").
%%
%% =============================================
%% IMPORTANT NOTICE:
%% See the copyright and distribution conditions below.
%%
%% Copyright 2011:
%%  Matthew R Kuhn
%%  Donald P. Shiley School of Engineering
%%  University of Portland
%%  5000 N Willamette Blvd
%%  Portland OR  97203
%%  USA
%%  E-mail: kuhn@up.edu
%%
%% This program can be redistributed and/or modified under the terms
%% of the LaTeX Project Public License Distributed, available from the CTAN
%% archives; either version 1.1 of the License, or any later version.
%% Most importantly, you may modify this file, but only if you rename it.
%%
%% This is a contributed file to the LaTeX2e system.
%%
%% This program consists of the files ascelike.cls and ascelike.bst
%%
%% -----------------------------------
%%
%% This document class produces bibliographic entries that roughly comply
%% with the guidelines of the American Society of Civil Engineers.
%% However, it was not produced by ASCE, its agents, or employees; nor 
%% is it in any way sanctioned or approved by that organization.
%%
%% Please contact Matthew Kuhn by email for suggestions, etc.
%%
%% -----------------------------------
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
    url
    volume
    year
  }
  {}
  { label.year label extra.label sort.label sort.year }

INTEGERS { output.state before.all mid.sentence after.sentence after.block 
           after.quote after.quote.comma 
           t.state o.state title.already org.already }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
  #4 'after.quote :=
  #5 'after.quote.comma :=
  #0 't.state :=
  #0 'o.state :=
  #1 'title.already :=
  #2 'org.already :=
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
            { output.state after.quote =
                { ".''\ " * write$ }
                { output.state after.quote.comma =
                    { ",''\ " * write$ }
	            { add.period$ " " * write$ }
                  if$
                }
              if$
            }
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

%					apalike needs this function because
%					the year has special punctuation;
%					apalike ignores the month
FUNCTION {output.year.check}
{ year empty$
    { "empty year in " cite$ * warning$ }
    { write$
      " (" year * extra.label * ")" *
      after.sentence 'output.state :=
    }
  if$
}

FUNCTION {fin.entry}
{ add.period$
  write$
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

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
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

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{vv~}{ll}{, jj}{, f.}" format.name$ 't :=   % last name first
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
	    { numnames #2 >
		{ "," * }
		'skip$
	      if$
	      t "others" =
		{ " et~al.\@" * }
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

FUNCTION {my.full.label}
{ 's :=
  #1 'nameptr :=               % nameptr = 1;
  s num.names$ 'numnames :=    % numnames = num.name$(s);
  numnames 'namesleft :=
    { namesleft #0 > }

    { s nameptr "{vv~}{ll}" format.name$ 't :=  % get the next name
      nameptr #1 >
        { namesleft #1 >
              { ", " * t * }
              { numnames #2 >
                  { "," * }
                  'skip$
                if$
                t "others" =
                    { " et~al.\@" * }
                    { " and " * t * } % from Chicago Manual of Style
                  if$
               }
               if$
             }
            't
        if$
        nameptr #1 + 'nameptr :=          % nameptr += 1;
        namesleft #1 - 'namesleft :=      % namesleft =- 1;
    }
  while$

}

FUNCTION {format.names.editors}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{f. }{vv}{ll}{, jj}" format.name$ 't :=   % last name first
      nameptr #1 >
        { namesleft #1 >
            { ", " * t * }
            { numnames #2 >
                { "," * }
                'skip$
              if$
              t "others" =
                { " et~al.\@" * }
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

FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.key}			% this function is just for apalike
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names.editors
      editor num.names$ #1 >
	{ ", eds." * }
	{ ", ed." * }
      if$
    }
  if$
}

INTEGERS { len.minus.1 }
STRINGS { u v }

FUNCTION {format.title}
{ title empty$
    { "" }
    { "``" title "t" change.case$ *
    }
%   { title 'u :=
%     u #1 #1 substring$ 'v :=
%     u text.length$ #1 - 'len.minus.1 :=
%     u #2 len.minus.1 substring$ 'u :=
%     u "l" change.case$ 'u :=
%     "``" v * u *
%   }
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

FUNCTION {format.btitle}
{ title emphasize
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

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "Vol." volume tie.or.space.connect
      series empty$
	'skip$
	{ " of " * series emphasize * }
      if$
      "Vol. and No." number either.or.check
    }
  if$
}

FUNCTION {format.howpublished}
{ howpublished empty$
  { "" }
  { "(" howpublished * ")" *
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
	{ "" pages n.dashify tie.or.space.connect }
	{ "" pages tie.or.space.connect }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume field.or.null
  number empty$
    'skip$
    { "(" number * ")" * *
      volume empty$
	{ "there's a number but no volume in " cite$ * warning$ }
	'skip$
      if$
    }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
	{ pop$ format.pages }
	{ ", " * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
	{ "Chapter" }
	{ type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
	'skip$
	{ ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
        { booktitle emphasize }
        { booktitle emphasize ", " * format.editors * }
%       { "In " booktitle emphasize * }
%      	{ "In " format.editors * ", " * booktitle emphasize * }
      if$
    }
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Report No." }
    { type }
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
  emphasize
}

FUNCTION {format.article.crossref}
{ "In"							% this is for apalike
  " \cite{" * crossref * "}" *
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
  "\cite{" * crossref * "}" *				% this is for apalike
}

FUNCTION {format.incoll.inproc.crossref}
{ "In"							% this is for apalike
  " \cite{" * crossref * "}" *
}

FUNCTION {format.lab.names}
{ 's :=
  s #1 "{vv~}{ll}" format.name$
  s num.names$ duplicate$
  #2 >
    { pop$ " et~al.\@" * }
    { #2 <
	'skip$
	{ s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " et~al.\@" * }
	    { " and " * s #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
	{ cite$ #1 #3 substring$ }
	'key					% apalike uses the whole key
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.organization.label}
%
% added - gnp. Provide label formatting by organization if editor is null.
%
{ editor empty$
    { organization empty$
        { key empty$
            { "no key, editor or organization in " cite$ * warning$
              cite$ #1 #3 substring$ }
            'key
          if$
        }
        { organization }
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
	{ key empty$
	    { cite$ #1 #3 substring$ }
	    'key				% apalike uses the whole key
	  if$
	}
	{ editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
	{ cite$ #1 #3 substring$ }
	'key			% apalike uses the whole key, no organization
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {author.key.organization.label}
%
% added - gnp. Provide label formatting by organization if author is null.
%
{ author empty$
    { key empty$
        { organization empty$
            { "no key, author or organization in " cite$ * warning$
              cite$ #1 #3 substring$ }
            'organization
          if$
        }
        { key }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {calc.label}
%
% Changed - GNP. See also author.organization.sort, editor.organization.sort
% Form label for BibTeX entry. The classification of which fields are used
% for which type of entry (book, inbook, etc.) are taken from alpha.bst.
% The change here from newapa is to also include organization as a
% citation label if author or editor is missing.
%
% Changed - MRK.  With ASCE's author-year citation style, there is no
% need for the first argument in the \citeauthoryear command.  The first
% argument in the chicago.bst style contains a full list of all
% author last names.  This argument is not used with ASCE's styel.
% We replace this argument with a null argument "{}",
% with the 2nd and 3rd arguments giving the citation author and the year.
% This change is made so that two separate articles, say by Huang, Bird, and
% Bell and by Haung, Bird, and Hendrich (both published in, say, 2009)
% will be given the citations Jones et al. (2009a) and Jones et al. (2009b).
% Without the change, both would be given the citations Jones et al. (2009).
% Thanks to Alec Gosse.
%
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
        'editor.key.organization.label
        { type$ "manual" =
            'author.key.organization.label
            'author.key.label
          if$
        }
      if$
    }
  if$

%------------- cut MRK -------------------
%  author empty$  % generate the full label citation information.
%    { editor empty$
%        { organization empty$
%            { type$ "manual" =
%                { title empty$
%                    { "no author, editor, organization or title in " 
%                      cite$ * warning$
%                      "??" 
%                    }
%                    { title }
%                  if$
%                }
%                { "no author, editor, or organization in " cite$ * warning$
%                  "??" }
%              if$
%            }
%            { organization }
%          if$
%        }
%        { editor my.full.label }
%      if$
%    }
%    { author my.full.label }
%  if$
%
%% leave label on the stack, to be popped when required.
%
%  "}{" * swap$ * "}{" *
%------------- paste MRK ---------------
   "}{"   swap$ * "}{" *
%------------- end MRK -----------------
%  year field.or.null purify$ #-1 #4 substring$ *
%
% save the year for sort processing afterwards (adding a, b, c, etc.)
%
  year field.or.null purify$ #-1 #4 substring$
  'label.year :=
}

FUNCTION {output.bibitem}
%{ newline$
%  "\bibitem[" write$
%  label write$
%  "]{" write$
%  cite$ write$
%  "}" write$
%  newline$
%  ""
%  before.all 'output.state :=
%}
{ newline$

  "\bibitem[\protect\citeauthoryear{" write$
  calc.label write$
  sort.year write$
  "}]{" write$

  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check 
  after.quote 'output.state :=
  crossref missing$
    { journal emphasize "journal" output.check
      organization output
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note empty$
    { }
    { before.all 'output.state :=
      "\ " note * output
    }
  if$
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
    }
    { format.authors output.nonnull
      crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  output.year.check				% special for apalike
  new.block
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      new.sentence
      format.howpublished output
      publisher "publisher" output.check
      address output
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  url empty$
    { }
    { mid.sentence 'output.state :=
      "$<$" url * "$>$" * output
%     format.url output
    }
  if$
  month empty$
    { }
    { before.all 'output.state :=
      "\ (" month * ")" * output
    }
  if$
  new.block
  note empty$
    { }
    { before.all 'output.state :=
      "\ " note * output
    }
  if$
% note output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check
  new.block
  howpublished output
  address output
  new.block
  note output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
    }
    { format.authors output.nonnull
      crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  output.year.check				% special for apalike
  new.block
  format.btitle "title" output.check
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      format.howpublished output
      publisher "publisher" output.check
      address output
      format.edition output
      format.chapter.pages output
    }
    { format.chapter.pages "chapter and pages" output.check
      new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  new.block
  note output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check
  after.quote 'output.state :=
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      format.howpublished output
      publisher "publisher" output.check
      address output
      format.edition output
      format.chapter.pages output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check
  after.quote 'output.state :=
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
      format.howpublished output
      organization output				% for apalike
      address output					% there's no year
%     new.sentence					% here so things
      publisher output					% are simpler
      format.pages output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
  url empty$
    { }
    { mid.sentence 'output.state :=
      "$<$" url * "$>$" * output
%     format.url output
    }
  if$
  month empty$
    { }
    { before.all 'output.state :=
      "\ (" month * ")" * output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  format.authors output
  #0 't.state :=
  #0 'o.state :=
  author empty$
    { format.editors output
      editor empty$
        { organization empty$
            { org.already 'o.state :=
              title empty$ 
                { title format.key output }
                { format.btitle "title" output.check
                  title.already 't.state :=
                }
              if$
            }
            { organization "organization" output.check 
              org.already 'o.state :=
            }
          if$
        }
        {  }
      if$
    }
    {  }
  if$
%   { format.editors "author and editor" output.check
%     editor format.key output
%   }
%   { format.authors output.nonnull
%     crossref missing$
%       { "author and editor" editor either.or.check }
%       'skip$
%     if$
%   }
% if$
% author format.key output				% special for
  output.year.check					% apalike
  new.block
  t.state title.already =
    {  }
    { format.btitle "title" output.check }
  if$
% format.btitle "title" output.check
  o.state org.already =
    { publisher address new.block.checkb
      publisher output
    }
    { organization address new.block.checkb
      organization output 
    }
  if$
% organization address new.block.checkb
% organization output
  address output
  format.edition output
  url empty$
    { }
    { mid.sentence 'output.state :=
      "$<$" url * "$>$" * output
%     format.url output
    }
  if$
  month empty$
    { }
    { before.all 'output.state :=
      "\ (" month * ")" * output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  new.block
  format.title "title" output.check
  after.quote 'output.state :=
  "M.S. thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output.nonnull
  url empty$
    { }
    { mid.sentence 'output.state :=
      "$<$" url * "$>$" * output
%     format.url output
    }
  if$
  month empty$
    { }
    { before.all 'output.state :=
      "\ (" month * ")" * output
    }
  if$
  new.block
  note output
  fin.entry
}
%{ output.bibitem
%  format.authors "author" output.check
%  author format.key output				% special for
%  output.year.check					% apalike
%  new.block
%  format.title "title" output.check
%  new.block
%  "Master's thesis" format.thesis.type output.nonnull
%  school "school" output.check
%  address output
%  new.block
%  note output
%  fin.entry
%}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check
  after.quote 'output.state :=
  journal emphasize "journal" output.check
  organization output
  url empty$
    { }
    { mid.sentence 'output.state := 
      "$<$" url * "$>$" * output 
%     format.url output
    }
  if$
  month empty$
    { }
    { before.all 'output.state :=
      "\ (" month * ")" * output 
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {format.url}
{ url empty$
    { "" }
    { "$<$" url * "$>$" *
    }
  if$
}

%{ output.bibitem
%  format.authors "author" output.check
%  author format.key output                              % special for
%  output.year.check                                     % apalike
%  new.block
%  format.title "title" output.check
%  after.quote 'output.state :=
%  crossref missing$
%    { journal emphasize "journal" output.check
%      organization output
%      format.vol.num.pages output
%    }
%    { format.article.crossref output.nonnull
%      format.pages output
%    }
%  if$
%  new.block
%  note output
%  fin.entry
%}




FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check
  after.quote 'output.state :=
  "Ph.D. thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output.nonnull
  url empty$
    { }
    { mid.sentence 'output.state :=
      "$<$" url * "$>$" * output 
%     format.url output
    }
  if$
  month empty$
    { }
    { before.all 'output.state :=
      "\ (" month * ")" * output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  format.editors output
  editor format.key output				% special for
  output.year.check					% apalike
  new.block
  format.btitle "title" output.check
  format.bvolume output
  format.number.series output
  address output				% for apalike
  new.sentence					% we always output
  organization output				% a nonempty organization
  publisher output				% here
  new.block
  note output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check
  after.quote 'output.state :=
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  url empty$
    { }
    { mid.sentence 'output.state :=
      "$<$" url * "$>$" * output
%     format.url output
    }
  if$
  month empty$
    { }
    { before.all 'output.state :=
      "\ (" month * ")" * output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% special for
  output.year.check					% apalike
  new.block
  format.title "title" output.check
  new.block
  note "note" output.check
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

%			There are three apalike cases: one person (Jones),
%			two (Jones and de~Bruijn), and more (Jones et~al.).
%			This function is much like format.crossref.editors.
%
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
      if$						% apalike uses initials
      s nameptr "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}" format.name$ 't := % <= here
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

FUNCTION {editor.sort}
{ editor empty$
    { key empty$
	{ "to sort, need editor or key in " cite$ * warning$
	  ""
	}
	{ key sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {author.organization.sort}
%
% added - GNP. Stack author or organization for sorting (from alpha.bst).
% Unlike alpha.bst, we need entire names, not abbreviations
%
{ author empty$
    { organization empty$
        { key empty$
            { "to sort, need author, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { organization sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.organization.sort}
%
% added - GNP. Stack editor or organization for sorting (from alpha.bst).
% Unlike alpha.bst, we need entire names, not abbreviations
%
{ editor empty$
    { organization empty$
        { key empty$
            { "to sort, need editor, organization, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { organization sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

%			apalike uses two sorting passes; the first one sets the
%			labels so that the `a's, `b's, etc. can be computed;
%			the second pass puts the references in "correct" order.
%			The presort function is for the first pass. It computes
%			label, sort.label, and title, and then concatenates.
%FUNCTION {presort}
%{ calc.label
%  label sortify
%  "    "
%  *
%  type$ "book" =
%  type$ "inbook" =
%  or
%    'author.editor.sort
%    { type$ "proceedings" =
%	'editor.sort
%	'author.sort
%      if$
%    }
%  if$
%  #1 entry.max$ substring$	% for
%  'sort.label :=		% apalike
%  sort.label			% style
%  *
%  "    "
%  *
%  title field.or.null
%  sort.format.title
%  *
%  #1 entry.max$ substring$
%  'sort.key$ :=
%}
FUNCTION {presort}
%
% Presort creates the bibentry's label via a call to calc.label, and then
% sorts the entries based on entry type. Chicago.bst adds support for
% including organizations as the sort key; the following is stolen from
% alpha.bst.
%
{ calc.label sortify % recalculate bibitem label
  year field.or.null purify$ #-1 #4 substring$ * % add year
  "    "
  *
  type$ "book" =
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
  #1 entry.max$ substring$        % added for newapa
  'sort.label :=                  % added for newapa
  sort.label                      % added for newapa
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

SORT		% by label, sort.label, title---for final label calculation

STRINGS { last.label next.extra }	% apalike labels are only for the text;

INTEGERS { last.extra.num }		% there are none in the bibliography

FUNCTION {initialize.extra.label.stuff}
{ #0 int.to.chr$ 'last.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
%
% Pass through all entries, comparing current entry to last one.
% Need to concatenate year to the stack (done by calc.label) to determine
% if two entries are the same (see presort)
%
{ last.label
  calc.label year field.or.null purify$ #-1 #4 substring$ * % add year
  #1 entry.max$ substring$ =     % are they equal?
     { last.extra.num #1 + 'last.extra.num :=
       last.extra.num int.to.chr$ 'extra.label :=
     }
     { "a" chr.to.int$ 'last.extra.num :=
       "" 'extra.label :=
       calc.label year field.or.null purify$ #-1 #4 substring$ * % add year
       #1 entry.max$ substring$ 'last.label := % assign to last.label
     }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
     'skip$
  if$
  label.year extra.label * 'sort.year :=
  extra.label 'next.extra :=
}

%FUNCTION {initialize.extra.label.stuff}	% and hence there is no `longest.label'
%{ #0 int.to.chr$ 'last.label :=
%  "" 'next.extra :=
%  #0 'last.extra.num :=
%}
%
%FUNCTION {forward.pass}
%{ last.label label =
%    { last.extra.num #1 + 'last.extra.num :=
%      last.extra.num int.to.chr$ 'extra.label :=
%    }
%    { "a" chr.to.int$ 'last.extra.num :=
%      "" 'extra.label :=
%      label 'last.label :=
%    }
%  if$
%}
%
%FUNCTION {reverse.pass}
%{ next.extra "b" =
%%    { "a" 'extra.label := }
%    'skip$
%  if$
%%  label extra.label * 'label :=
%  label.year extra.label * 'sort.year :=
%  extra.label 'next.extra :=
%}

EXECUTE {initialize.extra.label.stuff}

ITERATE {forward.pass}

REVERSE {reverse.pass}

%				Now that the label is right we sort for real,
%				on sort.label then year then title.  This is
%				for the second sorting pass.
FUNCTION {bib.sort.order}
{ sort.label
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

ITERATE {bib.sort.order}

SORT		% by sort.label, year, title---giving final bibliography order

FUNCTION {begin.bib}
{ preamble$ empty$				% no \etalchar in apalike
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{}" write$ newline$		% no labels in apalike
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
