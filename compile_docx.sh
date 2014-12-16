MDFILE=$1
DOCXFILE=${MDFILE//.md/.docx}

pandoc $MDFILE --to docx --from markdown+tex_math_single_backslash-implicit_figures --output $DOCXFILE --highlight-style tango --reference-docx template.docx
