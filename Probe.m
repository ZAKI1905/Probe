(* ::Package:: *)

(* :Title: Probe.m *)
(* :Author: M. Zakeri *)

BeginPackage["Probe`"]

Print["- \"Probe\" Package Loaded!\n- By M. Zakeri\n- Last updated on 11/24/2018.\n- Check \"Manual.nb\" for instructions."]

PlotLogGrid::"usage" = "PlotLogGrid[function -> F, grid-> {{X1min, X1max, X1res}, {X2min, X2max, X2res}, ...,{XNmin, XNmax, XNres}}, component-> i, contour-> c, reportcomponents-> {i_1, ..., i_m}, header -> {h_1, ..., h_n}, smooth -> True/False, valoutput -> \"Value_output.txt\",output -> \"output.txt\"]\n returns a list of the ND coordinates (contours)
along with the values for components i_1 to i_m of the function \"F\" on the grid specified by {X1,..., XN} for which the i-th component of F matches the given contour(s), i.e. F_i = c.";



Begin["`Private`"]


GetPoints[aplot_] := Cases[Normal[aplot], Line[pts_, ___] :> pts, Infinity]
 
 
 
SmoothPoints[points_, avg_] := Transpose[
     {MovingAverage[(#1[[1]] & ) /@ points, avg], 
      MovingAverage[(#1[[2]] & ) /@ points, avg],
      MovingAverage[(#1[[3]] & ) /@ points, avg],
      MovingAverage[(#1[[4]] & ) /@ points, avg]}]
 
 
 
exportTxt[TableOutput_, outfile_, OptionsPattern[format-> "C"]] := Module[{tbl1,tbl2}, 

 If[OptionValue[format]=="C", 
     tbl1=ToString[TableForm[Map[CForm, TableOutput,{2}]]];
	 tbl2=StringReplace[tbl1,"\n\n"->"\n"];, tbl2= TableOutput];
	
	
 Export[outfile, tbl2]; 
 Print[" ==> Table saved in:  '", outfile, "' ."]; ]; 
         
Grid2DAvg[list2D_, X_, Y_]:=Module[{base, \[Epsilon], avgpt},
base={Round[X],Round[Y]} ;
\[Epsilon]={X,Y}-base;

avgpt=(1-Abs[\[Epsilon][[1]]]-Abs[\[Epsilon][[2]]])list2D[[##]]&@@base 
      + Abs[\[Epsilon][[1]]]list2D[[base[[1]]+Sign[\[Epsilon][[1]]],base[[2]] ]] 
      +Abs[\[Epsilon][[2]]]list2D[[base[[1]],base[[2]]+Sign[\[Epsilon][[2]]] ]];
      
Return[avgpt];
]


(* ::Text:: *)
(* Finds the values of Func on the grid :*)


LogNGrid[Func_, grid_, OptionsPattern[LogNGrid]] := Module[{logiMin={}, logiMax={}, dlogi={}, Gtable, TableOutput,
cordVars, taskNums, cordVars10, tableVars, funcSimplify, ExcReg, ExcVal=0, ExcRegVar},
     
     logiMin= N[Log[10, grid[[;;,1]]]];
     logiMax= N[Log[10, grid[[;;,2]]]];
     dlogi = (logiMax - logiMin)/grid[[;;,3]];
      
       
      ExcRegVar = OptionValue[excludedregion]; 
      If[
          ExcRegVar === Null, ExcReg=EmptyRegion[Length[logiMin]],
          If[$VersionNumber<10, ExcReg=EmptyRegion[Length[logiMin]]; Print[" ==> Warning: \"RegionMember\" was added in version 10.\n ==> Exclusion region will be ignored."],
          If[Depth[ExcRegVar]==1, Print[" ==> Warning: No exclusion value provided, \"0\" will be used."],
           ExcReg=ExcRegVar[[1]]; ExcVal=ExcRegVar[[2]];];]
           
         ];
      
      ParallelEvaluate[locCounter=0;];
      progCounter = 0;
      SetSharedVariable[progCounter];
      taskNums = Times@@(#1+1&/@grid[[;;,3]]);
      cordVars=Table[ToExpression["logX"<>ToString[i]],{i,1, Length[logiMin]}];
      cordVars10=Table[ToExpression["10^logX"<>ToString[i]],{i,1, Length[logiMin]}];
      tableVars=Evaluate[Transpose[{cordVars, logiMin, logiMax, dlogi}]];
      
     
            
      funcSimplify[var_, func_, excreg_, excval_]:= If[$VersionNumber<10, func@@var, If[RegionMember[excreg,var], excval, func@@var]];    
        
      GTable[function_] := Monitor[
      
       ParallelTable[locCounter++; If[locCounter> taskNums/(2*$KernelCount), progCounter+=locCounter; locCounter=0;];
       
            funcSimplify[cordVars10, function, ExcReg, ExcVal], Evaluate[Sequence@@tableVars], DistributedContexts -> Automatic],
          
            Grid[ {{"Progress: ", ProgressIndicator[progCounter,{0, taskNums}], 
            "%"<>ToString[100*Round[progCounter/taskNums,0.001]]}} ,
              Frame -> All, Alignment -> {{Left, Center, Right}},ItemSize -> {{Automatic, Automatic, 3}}]
            ]; 
            
        
              
      TableOutput = AbsoluteTiming[GTable[Func]]; 
      
      Print[" ==> Time of computation: ", Round[First[TableOutput], 0.01], 
       " seconds!"];
         
      If[OptionValue[output] === Null, Null, exportTxt[ Prepend[Last[TableOutput], grid], OptionValue[output], format-> "Normal"]  ]; 
      
      Return[Last[TableOutput]]; 
      ]
 
Options[LogNGrid] = {output -> Null, excludedregion-> Null};


(* ::Text:: *)
(*Finds the plot points corresponding to the contour  Func[] == cont. :*)


(* ::Input::Initialization:: *)
PlotLogNGridFromTable::output = "The number of output files(`1`) should be the same as the contours(`2`)!";

PlotLogNGridFromTable[table_, cont_, grid_, OptionsPattern[{smooth -> False, output -> Null, component-> Null, reportcomponents-> Null, header-> Null}]] := Module[

     {logiMin={}, logiMax={}, dlogi={}, itoXi, GetContourPoints2D, 
      GridPlot2D, GridPlotSmooth2D, GPlot2D, GPlot, PltPoints},
      
     logiMin= N[Log[10, grid[[;;,1]]]];
     logiMax= N[Log[10, grid[[;;,2]]]];
     dlogi = (logiMax - logiMin)/grid[[;;,3]];
      
      itoXi[n_,ii_] := 10^(logiMin[[n]] + (ii - 1)*dlogi[[n]] );
      
      GetContourPoints2D[data_,contour_,itox_, jtoy_] := Module[{tempplot}, 
        tempplot = ListContourPlot[Transpose[data], Contours -> {{contour}}, 
           ContourShading -> None]; 
           ({itox[#1[[1]]], jtoy[#1[[2]]],#1[[1]],#1[[2]]} & ) /@ 
          Flatten[GetPoints[tempplot], 1]]; 
          
          GridPlot2D[points_, contour_, itox_, jtoy_] := GetContourPoints2D[points, contour, itox, jtoy]; 
       
      GridPlotSmooth2D[points_, contour_,itox_, jtoy_] := SmoothPoints[ GetContourPoints2D[points, contour, itox, jtoy], 3]; 
        
      GPlot2D[points_, contour_,itox_, jtoy_] := If[OptionValue[smooth] == True, 
        GridPlotSmooth2D[points, contour, itox, jtoy], GridPlot2D[points, contour, itox, jtoy]]; 
        

      GPlot[NDpoints_, contour_]:= Module[
                                                                                                  {VecSection2D, section2D,pltpts, finalpltpts={}, GridDim,DimVars, 
itoXiTAble,VecFuncVals,RepVecFuncVals,conturcoords,compRep, compRepSet}, 

GridDim=Length[logiMin];

(*DimVars=Table[{Symbol["m"<>ToString[i]],1, Dimensions[NDpoints][[i]]}, {i,1,GridDim-2}];*)

(* To generate: { {m1, 1, Res1+1}, {m2, 1, Res2+1},..., {m{GridDim-2}, 1, Res{GridDim-2}+1} }*)
DimVars=Table[{Symbol["ZakiM\[Alpha]"<>ToString[i]],1, grid[[;;,3]][[i]]+1}, {i,1,GridDim-2}];

(* { itoXi[1,m1],itoXi[2,m2],itoXi[3, m3],...,itoXi[GridDim-2, m{GridDim-2}] }*)
itoXiTAble= Table[itoXi[n,DimVars[[n,1]]],{n,1,GridDim-2}];

finalpltpts=ParallelTable[

	(* Vector 2D section of the ND values *)
     VecSection2D = Part[NDpoints, Evaluate[Sequence@@(#[[1]]&/@DimVars)]];

	(* Takes the contComp of the vector field, i.e. returns a scalar field for GPlot2D *)
	
section2D =If[OptionValue[component]=== Null,VecSection2D, #[[;;,OptionValue[component]]]&/@VecSection2D];


pltpts=GPlot2D[section2D, contour, itoXi[GridDim-1,#]&, itoXi[GridDim,#]&];
conturcoords={Evaluate[Sequence@@itoXiTAble],#[[1]],#[[2]]}&/@pltpts;

            
compRep=OptionValue[reportcomponents];



	If[compRep === Null,  conturcoords, 

                 If[Depth[compRep]==1, compRepSet={compRep}, compRepSet=compRep];
                  (*Returns NDpoints\[LeftDoubleBracket]m1, m2, m3, ..., mD\[RightDoubleBracket];*)

	     (*VecFuncVals=Part[VecSection2D, Round[#\[LeftDoubleBracket]3\[RightDoubleBracket]], Round[#\[LeftDoubleBracket]4\[RightDoubleBracket]] ]&/@pltpts;*)
                    VecFuncVals=Grid2DAvg[VecSection2D,#[[3]],#[[4]]]&/@pltpts;

	     If[
       Depth[VecFuncVals]>2, RepVecFuncVals =Table[#[[compRepSet[[i]]]], {i, Length[compRepSet]}]&/@VecFuncVals;,
           RepVecFuncVals =VecFuncVals;
       ];

Table[Flatten[Transpose[{conturcoords, RepVecFuncVals}][[i]]],{i, 1, Length[RepVecFuncVals]}]],

Evaluate[Sequence@@DimVars], DistributedContexts -> Automatic
                            ];

Return[ Flatten[finalpltpts,GridDim-2] ];   
];
        
      PltPoints =  GPlot[table, cont];

	If[OptionValue[header]=== Null,Null, PrependTo[PltPoints, OptionValue[header]];];
	
      If[OptionValue[output] === Null, Null, exportTxt[ PltPoints, OptionValue[output]]  ]; 
    
      Return[ PltPoints ]; 
      ]



 PlotLogGrid[OptionsPattern[PlotLogGrid]]:= Module[

     {logGridTable, logGridPoints, contourSet, contOpt, outputOpt, outputSet, inputgrid, infile, inputTable, Func, ExcReg, gridDims}, 
     
     inputgrid = OptionValue[grid];
     Func = OptionValue[function];
     infile = OptionValue[inputfile];
     ExcReg = OptionValue[excludedregion];
     
     If[inputgrid===  Null && infile===Null,
       
       For[i=1; gridDims={}; inputgrid={}; logGridTable ={}, i<=Length[Func], i++,
           AppendTo[logGridTable, LogNGrid[Func[[i,1]], Func[[i,2]], excludedregion-> ExcReg]]; 
           AppendTo[inputgrid, Func[[i,2]]];
           AppendTo[gridDims, Length[ Func[[i,2]] ] ];
          ];
       
       logGridTable= Outer[OptionValue[outeroperation],##]&@@ Flatten[{logGridTable, gridDims},1];
       (*logGridTable = Flatten[logGridTable, 1];*)
       inputgrid = Flatten[inputgrid, 1];
        If[OptionValue[valoutput] === Null, Null, exportTxt[ Prepend[logGridTable, inputgrid], OptionValue[valoutput], format-> "Normal"]  ]; ,
       
         (* inputgrid \[NotEqual] Null *)
         If[infile===Null,  logGridTable = LogNGrid[Func, inputgrid, output-> OptionValue[valoutput], excludedregion-> ExcReg]; ,
            inputTable = Import[infile,"Data"]//ToExpression; inputgrid = inputTable[[1]];  logGridTable = inputTable[[2;;]];
           ];
       ];
   
   
     
    contOpt=OptionValue[contour];
    outputOpt=OptionValue[output];
    
    If[
     contOpt=== Null, Return[logGridTable],
    
     If[
         outputOpt===Null, Null, If[
                                     Length[contOpt]!= Length[outputOpt], 
                                     Message[PlotLogNGridFromTable::output, Length[outputOpt], Length[contOpt]]
                                    ]  
       ];
     
     If[Depth[contOpt]==1, contourSet={contOpt}, contourSet=contOpt];
     If[Depth[outputOpt]==1, outputSet=Table[outputOpt, {dummy,1, Length[contourSet]}], outputSet=outputOpt];
    
     logGridPoints=PlotLogNGridFromTable[logGridTable, #[[1]], inputgrid, 
     output-> #[[2]], 
     smooth-> OptionValue[smooth],
     component-> OptionValue[component],  
     reportcomponents-> OptionValue[reportcomponents], 
     header-> OptionValue[header] ] &/@ Transpose[{contourSet, outputSet}];
     
     Return[logGridPoints]
     ]
     
 ]
     
Options[PlotLogGrid] = {contour-> Null, smooth -> False, output -> Null,
 valoutput-> Null, component-> Null,  reportcomponents-> Null, 
 inputfile-> Null, header-> Null, grid-> Null, function-> Null,
 excludedregion-> Null, outeroperation-> Null};


End[]
EndPackage[]
