#**********************************************************************************************#
#                                            MS Biomarkers
#  Copyright (C) 2024-2025 Osamah Alrouwab <rawab@uoz.edu.ly> Faculty of Medicine Zintan-university-Libya
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
#************************************************************************************************#

# ---------- Export DOCX (combined + individual) ----------
docx_path_both <- file.path(out_dir, "Subsection4_Tables.docx")

doc <- officer::read_docx()
doc <- officer::body_add_par(doc, "Subsection 4 – Multivariable modeling and predictors of disability", style = "heading 1")
doc <- officer::body_add_par(doc, "Table 5. Predictors of EDSS severity", style = "heading 2")
doc <- flextable::body_add_flextable(doc, value = ft5)
doc <- officer::body_add_par(doc, "")
doc <- officer::body_add_par(doc, "Table 6. Partial R² contributions", style = "heading 2")
doc <- flextable::body_add_flextable(doc, value = ft6)

print(doc, target = docx_path_both)

# individual tables
print(officer::read_docx() |> flextable::body_add_flextable(ft5),
      target = file.path(out_dir, "Table5_Ordinal_Logistic_EDSS.docx"))
print(officer::read_docx() |> flextable::body_add_flextable(ft6),
      target = file.path(out_dir, "Table6_PartialR2_EDSS.docx"))
