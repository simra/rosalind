// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.
open System
open System.Collections.Generic
open System.IO
open System.Text.RegularExpressions

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)


type NwckTree =
    | Leaf of string
    | Anon
    | Internal of (string*NwckTree list)

let ex = fun s -> new Exception(s)

let rec parseTree (str:string) =
    if str="" then raise (ex "Empty string")
    else if str.[0]<>'(' then raise (ex "Expected open brace")
    else
        let findClose () =
            let mutable o=0
            let mutable i=1
            while (i<str.Length && (str.[i]<>')' || o>0)) do
                if str.[i]='(' then o<- o+1
                else if str.[i]=')' then o<-o-1
                i<-i+1
            i
        let findLbl i = // TODO: I think we should never see close brace here.
            if i+1>=str.Length || str.[i+1]=';' || str.[i+1]=')' || str.[i+1]=',' then ("",i+1)
            else
                let mutable j=i+2
                while j<str.Length && str.[j]<>';' && str.[j]<>',' && str.[j]<>')' do j<-j+1
                let start=i+1
                let len=j-i-1
                //eprintfn "%s %d %d %d" str str.Length start len                
                (str.Substring(start,len),j+1)
        let i = findClose ()
        let (l,j) = findLbl i
        Internal(l,parseList (str.Substring(1,i-1)))// :: parseTree (str.Substring(j))
and parseList str =
    let commaSep ()=
        Seq.unfold(
            fun (o,i,last) ->
                if i>str.Length then None
                else if i=str.Length then
                    if o=0 then
                        Some(Some(str.Substring(last+1)),(0,i+1,i))
                    else
                        raise (ex "Unexpected end of string")
                else if o=0 && str.[i]=',' then
                    Some(Some(str.Substring(last+1,i-last-1)),(0,i+1,i))
                else if str.[i]='(' then 
                    Some(None,(o+1,i+1,last))
                else if str.[i]=')' then
                    Some(None,(o-1,i+1,last))
                else Some(None,(o,i+1,last))
            ) (0,0,-1)
        |> Seq.choose id
    commaSep() 
    |> Seq.map (fun s ->
        if s="" then Anon
        else if s.[0]='(' then parseTree s
        else Leaf(s)
    )
    |> List.ofSeq

let adjacencies _nwk =
    let d = new Dictionary<string,string list>()
    let rec adj nwk : unit =        
        let checkAnon s =
            if s="" then 
                let lbl=sprintf "anon_%d" (d.Count+1)
                d.[lbl]<-List.empty
                lbl
            else 
                if not (d.ContainsKey(s)) then d.[s]<-List.empty
                s

        match nwk with
        | Leaf(s) -> ()
        | Anon -> ()
        | Internal(l,n) ->
            let l' = checkAnon l                
            n
            |> Seq.iter 
                (fun nn ->
                    match nn with
                    | Leaf(s) -> 
                        let s'=checkAnon s
                        d.[s']<-l'::d.[s']
                        d.[l']<-s'::d.[l']
                    | Anon ->
                        let s'= checkAnon ""
                        d.[s']<-l'::d.[s']
                        d.[l']<-s'::d.[l']
                    | Internal(s,sub) ->
                        let s'=checkAnon s
                        d.[s']<-l'::d.[s']
                        d.[l']<-s'::d.[l']
                        adj (Internal(s',sub))
                
                )   
    adj _nwk
    d
// inefficient?  Proper handling of dead ends?
let distance (adj:Dictionary<string,string list>) a b =
    let visited=new HashSet<string>()
    let safeValue=adj.Count*10;
    let rec dist x =
        //eprintfn "Visiting %s" x        
        visited.Add(x)|> ignore        
        if x=b then 0
        else 
           // screams out for a computation expression     
           adj.[x]
           |> Seq.map                 
                (fun n ->
                    if not (visited.Contains(n)) then
              //          eprintfn "Trying %s" n
                        1 + (dist n)
                    else
                        safeValue)
           |> Seq.min 
                     
    dist a    

getData "nwck"
|> splitNewline
|> Array.ofSeq
|> fun s -> 
        seq {
            for i in [0..2..s.Length-1] do
                yield (parseTree s.[i]|>adjacencies,s.[i+1].Split(' '))
            }
//|> Seq.take 1
|> Seq.map (fun (adj,toks) -> distance adj toks.[0] toks.[1])
|> Seq.map string
|> String.concat " "
|> printfn "%s"    

// CTBL
let toEdges t =
    t
    |>adjacencies
    |>Seq.map (fun kvp -> seq { for v in kvp.Value do yield (kvp.Key,v)})
    |>Seq.concat

let isAnon s = Regex.IsMatch(s,@"anon_\d+")

let getCharRow edges (rem1,rem2) =
    edges
    |> Seq.fold 
        (fun (componentLookup:Map<string,int>,sets:Map<int,string list>,nextIx) (v1,v2) ->
            if (rem1=v1 && rem2=v2) || (rem1=v2 && rem2=v1) then 
          //      eprintfn "ignoring %s %s" v1 v2
                (componentLookup,sets,nextIx)
            else                    
                if Map.containsKey v1 componentLookup then
                    let c1=componentLookup.[v1]
                    if Map.containsKey v2 componentLookup then
                        let c2=componentLookup.[v2]
                        if c1=c2 then (componentLookup,sets,nextIx)
                        else
                            if List.length sets.[c2]>List.length sets.[c1] then
                                let outLookup = Seq.fold (fun m vi -> Map.add vi c2 m) componentLookup sets.[c1]
                                let outSets = sets|>Map.add c2 (List.concat [sets.[c2];sets.[c1]]) |> Map.remove c1
                                (outLookup,outSets,nextIx)
                            else 
                                let outLookup = Seq.fold (fun m vi -> Map.add vi c1 m) componentLookup sets.[c2]
                                let outSets = sets|>Map.add c1 (List.concat [sets.[c1];sets.[c2]]) |> Map.remove c2
                                (outLookup,outSets,nextIx)
                    else
                        (Map.add v2 c1 componentLookup, Map.add c1 (v2::sets.[c1]) sets, nextIx)
                else
                    if Map.containsKey v2 componentLookup then
                        let c2=componentLookup.[v2]
                        (Map.add v1 c2 componentLookup, Map.add c2 (v1::sets.[c2]) sets, nextIx)
                    else
                        let newSet=nextIx
                        (Map.add v1 newSet componentLookup|>Map.add v2 newSet, Map.add newSet [v1;v2] sets, newSet+1)
                        ) (Map.empty,Map.empty,0)
    |> fun (componentLookup,sets,_) ->
      //  eprintfn "Lbls: %d Sets: %d" (Seq.length componentLookup) (Seq.length sets)
        let ix =
            sets
            |> Map.toArray
            |> Array.map (fun (i,_) -> i)
            |> fun a -> seq { for i in [0..a.Length-1] do yield (a.[i],i)}
            |> Map.ofSeq
        
        componentLookup
        |> Map.toSeq
        |> Seq.filter (fun (s,_) -> not (isAnon s))
        |> Seq.map (fun (s,v) -> (s,ix.[v]))
        |> Seq.sortBy (fun (s,v)->s)
        |> Seq.map (fun (s,v)->v)
        |> Seq.map string
        |> String.concat ""
        
        

// seems there should be a faster (iterative) way to walk through the edges and relabel vertices
let ctbl t =
    let edges=toEdges t
    //edges  |> Seq.iter (fun (v1,v2)-> printfn "%s %s" v1 v2)
    edges
    |> Seq.filter (fun (v1,v2) -> isAnon v1 && isAnon v2) // ignore trivial edges. Note if we don't then the CC code isn't quite right (leafs get dropped).
    |> Seq.map (getCharRow edges)
    |> Seq.groupBy id
    |> Seq.map (fun (k,v)->k)
    |> Seq.iter (printfn "%s")


getData "ctbl"
|> splitNewline
|> Seq.exactlyOne
|> parseTree
|> ctbl


// NKEW.
// So the question is do I refactor the newick code above or rewrite it here...
// Refactoring was relatively fast and painless. Ideally we'd adjust parseWtTree to handle
// unweighted or weighted trees.
type WtNwckTree =
    | WtLeaf of (float*string)
    | WtAnon of float
    | WtInternal of (float*string*WtNwckTree list)


let rec parseWtTree (str:string) =
    if str="" then raise (ex "Empty string")
    else if str.[0]<>'(' then raise (ex "Expected open brace")
    else
        let findClose () =
            let mutable o=0
            let mutable i=1
            while (i<str.Length && (str.[i]<>')' || o>0)) do
                if str.[i]='(' then o<- o+1
                else if str.[i]=')' then o<-o-1
                i<-i+1
            i
        let getWt (s:string) =
            s.Split(':')
            |> fun toks -> 
                if toks.Length<2 then
                    eprintfn "Warning: no weight for string %s" s
                    0.,s
                else
                    (float toks.[1],toks.[0])
        let findLbl i = // TODO: I think we should never see close brace here.
            if i+1>=str.Length || str.[i+1]=';' || str.[i+1]=')' || str.[i+1]=',' then ("",i+1)
            else
                let mutable j=i+2
                while j<str.Length && str.[j]<>';' && str.[j]<>',' && str.[j]<>')' do j<-j+1
                let start=i+1
                let len=j-i-1
                //eprintfn "%s %d %d %d" str str.Length start len                
                (str.Substring(start,len),j+1)
        let i = findClose ()
        let (l,j) = findLbl i
        let (wt,lbl)=getWt l
        WtInternal(wt,lbl,parseWtList (str.Substring(1,i-1)))// :: parseTree (str.Substring(j))
and parseWtList str =
    let commaSep ()=
        Seq.unfold(
            fun (o,i,last) ->
                if i>str.Length then None
                else if i=str.Length then
                    if o=0 then
                        Some(Some(str.Substring(last+1)),(0,i+1,i))
                    else
                        raise (ex "Unexpected end of string")
                else if o=0 && str.[i]=',' then
                    Some(Some(str.Substring(last+1,i-last-1)),(0,i+1,i))
                else if str.[i]='(' then 
                    Some(None,(o+1,i+1,last))
                else if str.[i]=')' then
                    Some(None,(o-1,i+1,last))
                else Some(None,(o,i+1,last))
            ) (0,0,-1)
        |> Seq.choose id
    commaSep() 
    |> Seq.map (fun s ->
        //if s="" then Anon
        if s.[0]=':' then WtAnon(float (s.Substring(1)))
        else if s.[0]='(' then parseWtTree s
        else WtLeaf(float (s.Split(':').[1]),s.Split(':').[0])
    )
    |> List.ofSeq

let wtAdjacencies _nwk =
    let d = new Dictionary<string,(float*string) list>()
    let rec adj nwk : unit =        
        let checkAnon s =
            if s="" then 
                let lbl=sprintf "anon_%d" (d.Count+1)
                d.[lbl]<-List.empty
                lbl
            else 
                if not (d.ContainsKey(s)) then d.[s]<-List.empty
                s

        match nwk with
        | WtLeaf(w,s) -> ()
        | WtAnon(w) -> ()
        | WtInternal(w,l,n) ->
            let l' = checkAnon l                
            n
            |> Seq.iter 
                (fun nn ->
                    match nn with
                    | WtLeaf(w,s) -> 
                        let s'=checkAnon s
                        d.[s']<-(w,l')::d.[s']
                        d.[l']<-(w,s')::d.[l']
                    | WtAnon(w) ->
                        let s'= checkAnon ""
                        d.[s']<-(w,l')::d.[s']
                        d.[l']<-(w,s')::d.[l']
                    | WtInternal(w,s,sub) ->
                        let s'=checkAnon s
                        d.[s']<-(w,l')::d.[s']
                        d.[l']<-(w,s')::d.[l']
                        adj (WtInternal(w,s',sub))
                
                )   
    adj _nwk
    d
// inefficient?  Proper handling of dead ends?
let wtDistance (adj:Dictionary<string,(float*string) list>) a b =
    let visited=new HashSet<string>()
    let safeValue=float adj.Count*1e6;
    let rec dist x =
        //eprintfn "Visiting %s" x        
        visited.Add(x)|> ignore        
        if x=b then 0.
        else 
           // screams out for a computation expression     
           adj.[x]
           |> Seq.map                 
                (fun (w,n) ->
                    if not (visited.Contains(n)) then
              //          eprintfn "Trying %s" n
                        w + (dist n)
                    else
                        safeValue)
           |> Seq.min 
                     
    dist a    

getData "nkew"
|> splitNewline
|> Array.ofSeq
|> fun s -> 
        seq {
            for i in [0..2..s.Length-1] do
                yield (parseWtTree s.[i]|>wtAdjacencies,s.[i+1].Split(' '))
            }
//|> Seq.take 1
|> Seq.map (fun (adj,toks) -> wtDistance adj toks.[0] toks.[1])
|> Seq.map string
|> String.concat " "
|> printfn "%s"    

//ghbp
// requires validation
// seems that we can't guarantee consistency.
// so maybe we need some ordering on splits in order to build the correct tree.
// we could enumerate all possible combinations and choose the first that validates- too expensive for >9 rows.
// http://evolution.berkeley.edu/evolibrary/article/0_0_0/phylogenetics_07
type Species = { Name:string; Characters:Map<int,int> }
let isSuperset s1 s2 = // is s2 a subset of s1?
    Seq.zip s1 s2
    |> Seq.fold (fun b (c1,c2)-> b&&(c1>=c2)) true
let sortByPrecedence (ctbl:string list) =
    let isEqual s1 s2 =
        Seq.zip s1 s2
        |> Seq.fold (fun b (c1,c2) -> b&&(c1=c2)) true
    
    

    let complement s =
        s|>Seq.map (fun c -> 1-c)

    let ctbl' = 
        ctbl
        |> List.map (fun s -> s|>Seq.map (fun c -> string c |> int))

    let s0=List.head ctbl'
    ctbl'
    |> List.map 
        (fun si -> 
            if isSuperset s0 si || isSuperset si s0 then si 
            else 
                let t=complement si
                if (not (isSuperset s0 t) && (not (isSuperset t s0))) then
                    raise (new Exception "Failed to complement")
                t)
    |> Array.ofList 
    |> fun a ->
        for i in [a.Length-1..-1..0] do
            for j in [0..i-1] do
                if isSuperset a.[j] a.[i] then
                    let tmp = a.[i]
                    a.[i]<-a.[j]
                    a.[j]<-tmp
        a

(*
let reOrder s1 s2 =
    if isSuperset s1 s2 then
        ('A',s1,s2)
    else if isSuperset s2 s1 then
        ('B',s2,s1)
    else if isSuperset s1 (complement s2) then
        ('C',s1,complement s2)
    else if (isSuperset (complement s2) s1) then
        ('D',complement s2,s1)
    else ('E',s1,s2)


    let compareCharacters (s1:int seq) (s2:int seq)=
        Seq.zip s1 s2
        |> Seq.map (fun (c1,c2) -> (if c1=1 then (c1-c2) else 0), (if c2=1 then c2-c1 else 0))
        |> Seq.reduce (fun (a,b) (c,d) -> (a+c,b+d))
        *)

// 
let isSuperset s1 s2 = // is s2 a subset of s1?
    Seq.zip s1 s2
    |> Seq.fold (fun b (c1,c2)-> b&&(c1>=c2)) true

let complement s =
    s|>Seq.map (fun c -> 1-c)
let complementStr s =
    s|>Seq.map (fun c -> if c='1' then "0" else "1")|> String.concat ""

let toInts s = s|>Seq.map (fun c -> string c |> int)

let makePhylogeny tree splits =
    match splits with
    [] -> tree
    head::tail ->
        match tree with
        | Anon -> Leaf(head)
        | Leaf(l2) ->
            let s1=toInts head
            let s2=toInts l2
            if isSuperset s1 s2 then
                Internal(head,[Leaf(l2);Leaf(l2)])
            else if isSuperset s2 s1 then
                Internal(l2,[Leaf(head);Leaf(head)])
            else if isSuperset s1 (complement s2) then
                Internal(head,[Leaf(complementStr l2);Leaf(complementStr l2)])
            else // must be superset s2c over s1
                Internal(complement l2,[Leaf(head);Leaf(head)])
        | Internal(s,tlist) ->
            
getData "ghbp" 
|> splitNewline
|> fun arr ->
    let taxa=arr.[0].Split(' ')
    arr.[1..]
    |> List.ofArray
    |> List.fold makePhylogeny Anon
    |> sortByPrecedence //Array.map (fun s -> s|>Seq.map (fun c -> string c |> int)) 
    |> Array.ofSeq
    |> fun a -> eprintfn "%A" (isSuperset a.[0] a.[1]); a
    |> Seq.map (fun s -> s|>Seq.map string|>String.concat "")
    |> Seq.iter (printfn "%s")
    
    
    |> fun a -> 
        seq {
            for i in [0..a.Length-1] do
                for j in [0..a.Length-1] do
                //printfn "%d %d %A %A %A" i j a.[i] a.[j] (compareCharacters a.[i] a.[j])
                reOrder a.[i] a.[j]
         }
         |>Seq.filter (fun (s,_,_) -> s='E') //printfn "%c" s//i j (isSubset a.[i] a.[j])
         |> printfn "%A"

let prepCtbl (taxa:string[]) (ctbl:string list) =
    ctbl
    |> Seq.mapi (fun i s -> s|>Seq.mapi (fun j c -> taxa.[j],(i,string c|>int)))
    |> Seq.concat
    |> Seq.groupBy (fun (t,(k,v))->t)
    |> Seq.map (fun (t,s)-> t,s|>Seq.map (fun (t,(k,v))->(k,v))|>Map.ofSeq)
    |> Seq.map (fun (t,m)-> {Name=t; Characters=m})
    |> List.ofSeq

let species=
    getData "ghbp" 
    |> splitNewline
    |> fun arr ->
        let taxa=arr.[0].Split(' ')
        arr.[1..]
        |> List.ofSeq
        |> prepCtbl taxa


let splitTaxa species =
    let splits =
        species
        |> Seq.map (fun s -> s.Characters)
        |> Seq.map (fun m -> m|>Map.toSeq)
        |> Seq.concat
        |> Seq.groupBy (fun (k,v)->k)
        |> Seq.map (fun (k,s)-> (k,s|>Seq.countBy (fun (k,v)->v)|>Map.ofSeq))
        |> Seq.filter (fun (k,m)-> (Map.containsKey 0 m && Map.containsKey 1 m))

    let splitKey =
        if (Seq.isEmpty splits) then
            -1
        else 
            splits
            |> Seq.maxBy (fun (k,m)-> if (Map.containsKey 1 m) then m.[1] else 0)
            |> fun (k,m)->k
    
    let dropCharacter k s =
        s
        |> Seq.map (fun sp ->  {sp with Characters=Map.remove k sp.Characters})
        
    if splitKey>=0 then 
        species
        |> Seq.groupBy (fun s -> 
                s.Characters.[splitKey])
        |> Seq.map (fun (c,s)-> s)
        |> Seq.map (dropCharacter splitKey)
    else 
        seq { yield species }

let printSpecies species =
    species
    |> Seq.map (fun s -> s.Name)
    |> String.concat ";"
            

let rec buildPhylo species =
    let result = 
        splitTaxa species
        |> Array.ofSeq
        
    if result.Length=1 then
        eprintfn "%s" (species |> printSpecies)
        result.[0]
        |> Seq.map (fun s -> s.Name)
        |> String.concat ";"
        |> Leaf
    else
        Internal("",[buildPhylo result.[0];buildPhylo result.[1]])
        
(*
let rec splitTree (taxa:string[]) (nwck:NwckTree) (splitLine:string) =
    let split =
        splitLine
        |> Seq.mapi (fun i c -> taxa.[i],c) 
        |> Seq.groupBy (fun (t,c)->c)
        |> Map.ofSeq
    let splitLookup =
        splitLine
        |> Seq.mapi (fun i c -> taxa.[i],c)
        |> Map.ofSeq
        
    match nwck with
    | Anon ->
        let s1=
            split.['0']
            |> Seq.map (fun (t,c)->t)
            |> String.concat ";"
        let s2=
            split.['1']
            |> Seq.map (fun (t,c)->t)
            |> String.concat ";"
        Internal("",[Leaf(s1);Leaf(s2)])
    | Leaf(s) ->
        let result =
            s.Split(';')
            |> Seq.groupBy (fun s -> splitLookup.[s])
            |> Seq.map (fun (c,s) -> s)
            |> Seq.map (String.concat ";")
            |> Array.ofSeq
        if result.Length>1 then
            Internal("",[Leaf(result.[0]);Leaf(result.[1])])
        else
            Leaf(result.[0])
    | Internal(s,children) ->        
        Internal(s,children|>Seq.map (fun c -> splitTree taxa c splitLine)|>List.ofSeq)
 *)
let rec finalizeTree (t:NwckTree) =
    match t with
    | Anon -> Anon
    | Leaf(s) -> 
        s.Split(';')
        |> fun a -> 
            if (a.Length=1) then Leaf(s)
            else Internal("",[Leaf(a.[0]);finalizeTree (Leaf(String.concat ";" a.[1..]))]) 
    | Internal(s,leaves) ->
        Internal(s,leaves|>List.map finalizeTree)    

let unrootTree (t:NwckTree) =
    // the topmost tree is an internal with no parent- degree two.
    // take one of the two branches and merge it with the other.
    match t with
    | Internal(s,leaves) ->
        leaves
        |> Array.ofList
        |> fun a ->
           // if (a.Length!=2) then eprintfn "FAIL"
            let t0 = a.[0]
            match a.[1] with
            | Internal (s',l2) -> Internal(s',t0::l2)
            | _ -> 
                eprintfn "FAILED TO unroot"
                t
    | _ -> 
        eprintfn "FAILED TO unroot2"
        t



let rec printTree t =
    match t with
    | Anon -> ""
    | Leaf(s) -> s
    | Internal(s,leaves) ->
        leaves
        |> Seq.map printTree
        |> String.concat ","
        |> fun l -> sprintf "(%s)%s" l s

let rec addSpecies ctbl (v,e) s =
    let mkNew vertices =
        let vCount = Seq.length vertices
        sprintf "V%d" vCount :: vertices

    match e with
    | [] -> "Empty tree."
    | (v1,v2) :: tail ->
        let (vnew,v')=mkNew v
        let t2 =
            (s,vnew) :: (v1, vnew) :: (v2,vnew) :: tail  
        if isValid t2 then
            t2
        else
            let t2 = head :: (addSpecies ctbl (v,tail) s)
            if isValid t2 then 
                t2
            else
                raise (new Exception("Can't find valid addition"))

getData "ghbp" 
|> splitNewline
|> fun arr ->
    let taxa=arr.[0].Split(' ')
    arr.[1..]
    |> List.ofSeq
    |> prepCtbl taxa
    |> buildPhylo     
|> finalizeTree
|> unrootTree
|> printTree
|> printfn "%s;"

[<EntryPoint>]
let main argv = 
    printfn "%A" argv
    0 // return an integer exit code
