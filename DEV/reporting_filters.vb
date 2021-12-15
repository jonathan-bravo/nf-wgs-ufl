' Filter SNPs
'
' This macro applies a filter to the SNPs based on ClinVar, CADD, REVEL, and homozygosity
'
Sub FilterSNPs()
    Cells(1,20).Value = "Filter"
    Dim X As Long
    NumRows = Range("A1", Range("A1").End(xldown)).Rows.Count
    Range("A1").Select
    For X = 2 To NumRows
        If Cells(X,15).Value = 0 And Cells(X,16).Value < 20 And Cells(X,17).Value <> "Pathogenic" And Cells(X,17) <> "Pathogenic/Likely_pathogenic" And Cells(X,12).Value <> "homozygous" Then
            Cells(X,21) = "DON'T KEEP"
        ElseIf Cells(X,16).Value < 10 And Cells(X,15).Value < 0.5 And Cells(X,17).Value <> "Pathogenic" And Cells(X,17) <> "Pathogenic/Likely_pathogenic" And Cells(X,12).Value <> "homozygous" Then
            Cells(X,21) = "DON'T KEEP"
        Else
            Cells(X,21) = "KEEP"
        End If
        Dim Y As Long
        Dim ReadsArray(4) As Integer
        OriginalData = Replace(Cells(X,13), ":", ",")
        CommaData = Split(OriginalData, ",")
        For Y = 0 To 3
            ReadsArray(Y) = Cint(CommaData(Y))
        Next Y
        TotalSupportingReads = Application.Sum(ReadsArray)
        ' If the ClinVar significance is Pathogenic, then the variant is highlighted red
        If Cells(X,17) = "Pathogenic" Or Cells(X,17) = "Pathogenic/Likely_pathogenic" Then
            Rows(X).Interior.Color = RGB(255,199,206)
            Rows(X).Font.Color = RGB(156,0,6)
        ' If the CADD score is >= 20 and REVEL score is >= 0.5 then the variant is highlighted yellow
        ElseIF Cells(X,16).Value >= 20 AND Cells(X,15).Value >= 0.5 And TotalSupportingReads >= 10 Then
            Rows(X).Interior.Color = RGB(255,255,153)
            Rows(X).Font.Color = RGB(147,90,30)
        ' If the variant is homozygous then it is highlighted green
        ElseIf Cells(X,12) = "homozygous" And TotalSupportingReads >= 10 Then
            Rows(X).Interior.Color = RGB(204,238,206)
            Rows(X).Font.Color = RGB(34,95,0)
        Else
            Rows(X).Interior.Color = RGB(255,255,255)
            Rows(X).Font.Color = RGB(0,0,0)
        End IF
        ActiveCell.Offset(1, 0).Select
    Next X
    Range("A2").AutoFilter _
        Field:=21, _
        Criteria1:="KEEP"
    Range("U:U").EntireColumn.Hidden = TRUE
End Sub

' Filter SVs
'
' This macro applies a filter to the SVs based on ClinVar, CADD, and homozygosity
'
Sub FilterSVs()
    Cells(1,19).Value = "Filter"
    Dim X As Long
    NumRows = Range("A1", Range("A1").End(xldown)).Rows.Count
    Range("A1").Select
    For X = 2 To NumRows
        If Cells(X,14).Value < 10 And Cells(X,16).Value <> "Pathogenic" And Cells(X,16) <> "Pathogenic/Likely_pathogenic" And Cells(X,12).Value <> "homozygous" Then
            Cells(X,20) = "DON'T KEEP"
        Else
            Cells(X,20) = "KEEP"
        End If
        Dim Y As Long
        Dim ReadsArray(4) As Integer
        OriginalData = Replace(Cells(X,13), ":", ",")
        CommaData = Split(OriginalData, ",")
        For Y = 0 To 3
            ReadsArray(Y) = Cint(CommaData(Y))
        Next Y
        TotalSupportingReads = Application.Sum(ReadsArray)
        If Cells(X,16) = "Pathogenic" Or Cells(X,16) = "Pathogenic/Likely_pathogenic"Then
            Rows(X).Interior.Color = RGB(255,199,206)
            Rows(X).Font.Color = RGB(156,0,6)
        ElseIf Cells(X,14).Value >= 20 And TotalSupportingReads >= 10 Then
            Rows(X).Interior.Color = RGB(255,255,153)
            Rows(X).Font.Color = RGB(147,90,30)
        ElseIf Cells(X,12) = "homozygous" And TotalSupportingReads >= 10 Then
            Rows(X).Interior.Color = RGB(204,238,206)
            Rows(X).Font.Color = RGB(34,95,0)
        Else
            Rows(X).Interior.Color = RGB(255,255,255)
            Rows(X).Font.Color = RGB(0,0,0)
        End IF
        ActiveCell.Offset(1, 0).Select
    Next X
    Range("A2").AutoFilter _
        Field:=20, _
        Criteria1:="KEEP"
    Range("T:T").EntireColumn.Hidden = TRUE
End Sub

' Filter EXPs
'
' This macro applies a filter to the expansions based on normal and affected ranges
'
Sub FilterEXPs()
    
End Sub