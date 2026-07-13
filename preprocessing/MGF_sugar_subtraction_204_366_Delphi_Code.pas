// Start Tn and T subtraction. Recognize 204 and 366 and separate data. 

procedure TTForm1.nandTSubtraction1Click(Sender: TObject);
var   F_Import, F_Temp, F_Export1, F_Export2, F_Export3, F_Export4, F_Export5, F_Export6, F_Export7, F_Export8, F_Export9, F_Export10, F_Export11, F_Export12, F_Export13: TextFile;
      Pos_start, Pos_end, NCol, i, k, Flag, Flag1, z, PrecPos, RestPos:integer;
      S, S1, S2, S3, S4, St1, St2, Word, S_Export:String;
      Words:TStringList;
      Precursor:real;



begin
 OpenDialog1.Filter:='.mgf|*.mgf';
 if not OpenDialog1.Execute then
 Exit;
 AssignFile(F_Import,OpenDialog1.FileName);
 St1:=OpenDialog1.FileName;
 St2:=St1;
 Delete(St2,Pos('.mgf',St1),4);
 AssignFile(F_Export1,St2+'_(1G)_1xTn.mgf');
 AssignFile(F_Export2,St2+'_(2G)_2xTn.mgf');
 AssignFile(F_Export3,St2+'_(3G)_3xTn.mgf');
 AssignFile(F_Export4,St2+'_(1G)_1xT.mgf');
 AssignFile(F_Export5,St2+'_(2G)_1xTn_1xT.mgf');
 AssignFile(F_Export6,St2+'_(3G)_2xTn_1xT.mgf');
 AssignFile(F_Export7,St2+'_(2G)_2xT.mgf');
 AssignFile(F_Export8,St2+'_(3G)_1xTn_2xT.mgf');
 AssignFile(F_Export9,St2+'_(3G)_3xT.mgf');
 AssignFile(F_Export10,St2+'_(4G)_4xTn.mgf');
 AssignFile(F_Export11,St2+'_(4G)_4xT.mgf');
 AssignFile(F_Export12,St2+'_(5G)_5xTn.mgf');
 AssignFile(F_Export13,St2+'_(5G)_5xT.mgf');

 Rewrite(F_Export1);
 Rewrite(F_Export2);
 Rewrite(F_Export3);
 Rewrite(F_Export4);
 Rewrite(F_Export5);
 Rewrite(F_Export6);
 Rewrite(F_Export7);
 Rewrite(F_Export8);
 Rewrite(F_Export9);
 Rewrite(F_Export10);
 Rewrite(F_Export11);
 Rewrite(F_Export12);
 Rewrite(F_Export13);
 try
 Reset (F_Import);

  except ShowMessage ('Can not open the file'+OpenDialog1.FileName);
 Exit;
 end;
  Words:=TStringList.Create;
  Words.LoadFromFile(OpenDialog1.FileName);
    Pos_start:=0;
  Pos_end:=0;
  i:=0;
  k:=0;
  Flag:=0;
  Flag1:=0;

 ProgressBar1.Min:=0;
 ProgressBar1.Max:=Words.Count-1;

 WriteLn(F_Export1,'MASS=Monoisotopic');
 WriteLn(F_Export2,'MASS=Monoisotopic');
 WriteLn(F_Export3,'MASS=Monoisotopic');
 WriteLn(F_Export4,'MASS=Monoisotopic');
 WriteLn(F_Export5,'MASS=Monoisotopic');
 WriteLn(F_Export6,'MASS=Monoisotopic');
 WriteLn(F_Export7,'MASS=Monoisotopic');
 WriteLn(F_Export8,'MASS=Monoisotopic');
 WriteLn(F_Export9,'MASS=Monoisotopic');
 WriteLn(F_Export10,'MASS=Monoisotopic');
 WriteLn(F_Export11,'MASS=Monoisotopic');
 WriteLn(F_Export12,'MASS=Monoisotopic');
 WriteLn(F_Export13,'MASS=Monoisotopic');
 


 For i:=0 to Words.Count-1 do
  begin
      S:=Words.Strings[i];
      ProgressBar1.Position:=i;
      if Pos('BEGIN',S)>0 Then Pos_start:=i;
      if Pos('204.08',S)=1 Then Flag:=1;
      if Pos('366.13',S)=1 Then Flag1:=1;        // Toggle to recognize 366
      if Pos('366.14',S)=1 Then Flag1:=1;	// Toggle to recognize 366
      if Pos('END',S)>0 Then Pos_end:=i;

      if Pos('PEPMASS',S)>0 Then
                begin
                  PrecPos:=Pos('=',S)+1;
                  S1:=Copy(S,1,PrecPos-1);
                  S2:=S;
                  Delete(S2,1,PrecPos-1);
                  RestPos:=Pos(' ',S2)+1;
                  S3:=S2;
                  Delete(S3,1,RestPos-2); //Intensisty
                  S2:=Copy(S2,1,RestPos-2); //Precursor ion
                  Precursor:=StrToFloat(S2);
                end;
       if Pos('CHARGE',S)>0 Then
                begin
                  S4:=Copy(S,1,Length(S)-1);
                  Delete(S4,1,Pos('=',S));
                  z:=StrToInt(S4);          //Charge
                end;

      If (Pos_start>0) AND (Flag=1) AND (Flag1=0) AND (Pos_end>0) Then
        begin
          for k:=Pos_start to Pos_end do
            begin
              S:=Words.Strings[k];
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-203.07937)>=300) Then WriteLn(F_Export1,S1+FloatToStr(Precursor-203.07937/z)+S3); //Tn
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-203.07937)>=300) Then WriteLn(F_Export1,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-406.15874)>=300) Then WriteLn(F_Export2,S1+FloatToStr(Precursor-406.15874/z)+S3); //2Tn
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-406.15874)>=300) Then WriteLn(F_Export2,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-609.23811)>=300) Then WriteLn(F_Export3,S1+FloatToStr(Precursor-609.23811/z)+S3); //3Tn
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-609.23811)>=300) Then WriteLn(F_Export3,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-812.31748)>=300) Then WriteLn(F_Export10,S1+FloatToStr(Precursor-812.31748/z)+S3); //4Tn
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-812.31748)>=300) Then WriteLn(F_Export10,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-1015.39685)>=300) Then WriteLn(F_Export12,S1+FloatToStr(Precursor-1015.39685/z)+S3); //5Tn
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-1015.39685)>=300) Then WriteLn(F_Export12,S);
            end;
        end;
      If (Pos_start>0) AND (Flag=1) AND (Flag1=1) AND (Pos_end>0) Then       //Flag1 =0 ignore 366. Flag1=1 accept 366
        begin
          for k:=Pos_start to Pos_end do
            begin
              S:=Words.Strings[k];

              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-365.13219)>=300) Then WriteLn(F_Export4,S1+FloatToStr(Precursor-365.13219/z)+S3);   //T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-365.13219)>=300) Then WriteLn(F_Export4,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-568.21157)>=300) Then WriteLn(F_Export5,S1+FloatToStr(Precursor-568.21157/z)+S3); //Tn_T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-568.21157)>=300) Then WriteLn(F_Export5,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-771.29094)>=300) Then WriteLn(F_Export6,S1+FloatToStr(Precursor-771.29094/z)+S3); //2Tn_T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-771.29094)>=300) Then WriteLn(F_Export6,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-730.2644)>=300) Then WriteLn(F_Export7,S1+FloatToStr(Precursor-730.2644/z)+S3);   //2T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-730.2644)>=300) Then WriteLn(F_Export7,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-933.34377)>=300) Then WriteLn(F_Export8,S1+FloatToStr(Precursor-933.34377/z)+S3); //Tn_2T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-933.34377)>=300) Then WriteLn(F_Export8,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-1095.3966)>=300) Then WriteLn(F_Export9,S1+FloatToStr(Precursor-1095.3966/z)+S3); //3T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-1095.3966)>=300) Then WriteLn(F_Export9,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-1460.5288)>=300) Then WriteLn(F_Export11,S1+FloatToStr(Precursor-1460.5288/z)+S3); //4T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-1460.5288)>=300) Then WriteLn(F_Export11,S);
              if (Pos('PEPMASS',S)>0) AND ((Precursor*z-1825.66095)>=300) Then WriteLn(F_Export13,S1+FloatToStr(Precursor-1825.66095/z)+S3); //5T
              if (Pos('PEPMASS',S)=0) AND ((Precursor*z-1825.66095)>=300) Then WriteLn(F_Export13,S);
            end;
        end;


     If (Pos_start>0) And (Flag=1) AND (Pos_end>0) Then
        begin
          Pos_start:=0;
          Pos_end:=0;
          Flag:=0;
          Flag1:=0;
        end;

     If (Pos_start>0) And (Flag=0) AND (Pos_end>0) Then
        begin
          Pos_start:=0;
          Pos_end:=0;
          Flag:=0;
          Flag1:=0;
       end;
    If (Pos_start>0) And (Flag1=1) AND (Pos_end>0) Then
        begin
          Pos_start:=0;
          Pos_end:=0;
          Flag:=0;
          Flag1:=0;
        end;

     If (Pos_start>0) And (Flag1=0) AND (Pos_end>0) Then
        begin
          Pos_start:=0;
          Pos_end:=0;
          Flag:=0;
          Flag1:=0;
       end;

 end;
 
 CloseFile(F_Import);
 CloseFile(F_Export1);
 CloseFile(F_Export2);
 CloseFile(F_Export3);
 CloseFile(F_Export4);
 CloseFile(F_Export5);
 CloseFile(F_Export6);
 CloseFile(F_Export7);
 CloseFile(F_Export8);
 CloseFile(F_Export9);
 CloseFile(F_Export10);
 CloseFile(F_Export11);
 CloseFile(F_Export12);
 CloseFile(F_Export13);
 ProgressBar1.Position:=0;
end;
// End Tn and T subtraction. Recognize 204 and 366 and separate data
