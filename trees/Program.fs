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

let rec C n k =
    if k=0L || k=n then 1L
    else if k>=(n/2L+1L) then C n (n-k)
    else 
        [1L..k]
        |> Seq.map (fun i -> (float (n-k+i))/(float i))
        |> Seq.reduce (*)
        |> int64


open System.Collections.Generic
 // some handy dynamic programming snippets from http://www.fssnip.net/8P
 /// The function creates a function that calls the argument 'f'
 /// only once and stores the result in a mutable dictionary (cache)
 /// Repeated calls to the resulting function return cached values.
let memoize f =    
   // Create (mutable) cache that is used for storing results of 
   // for function arguments that were already calculated.
    let cache = new Dictionary<_, _>()
    (fun x ->
       // The returned function first performs a cache lookup
       let succ, v = cache.TryGetValue(x)
       if succ then v else 
         // If value was not found, calculate & cache it
         let v = f(x) 
         cache.Add(x, v)
         v)


type NwckTree =
    | Leaf of string
    | Anon
    | Internal of (string*NwckTree list)

let ex = fun s -> new Exception(s)

let rec printTree t =
    match t with
    | Anon -> ""
    | Leaf(s) -> s
    | Internal(s,leaves) ->
        leaves
        |> Seq.map printTree
        |> String.concat ","
        |> fun l -> sprintf "(%s)%s" l s


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
    


getData "ctbl"
|> splitNewline
|> Seq.exactlyOne
|> parseTree
|> ctbl
|> Seq.iter (printfn "%s")

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


// cntq
// could this just be C(n,4)?
// Trick question!
getData "cntq"
|> splitNewline
|> Array.ofSeq
|> fun a -> 
    let n = int64 a.[0]
    n*(n-1L)*(n-2L)*(n-3L)/(24L)
|> fun x -> x%1000000L
|> printfn "%d"

// eubt
// enumerate all the ways to split edges and insert the leaf.
// messy but effective.
let rec enumerateAddToTree t (leaf:string) =
    seq {
        match t with 
        | Internal(s,l) ->
            let a = l|>Array.ofList
            // can we make this more succinct?            
            if a.Length=2 then
                yield Some(Internal(s,[Internal("",[a.[0];Leaf(leaf)]);a.[1]]))
                yield Some(Internal(s,[a.[0];Internal("",[a.[1];Leaf(leaf)])]))
                yield! 
                    enumerateAddToTree a.[0] leaf
                    |> Seq.choose id
                    |> Seq.map (fun t' -> Some(Internal(s,[t';a.[1]])))
                yield!
                    enumerateAddToTree a.[1] leaf
                    |> Seq.choose id
                    |> Seq.map (fun t' -> Some(Internal(s,[a.[0];t'])))
            else
                yield Some(Internal(s,[Internal("",[a.[0];Leaf(leaf)]);a.[1];a.[2]]))
                yield Some(Internal(s,[a.[0];Internal("",[a.[1];Leaf(leaf)]);a.[2]]))
                yield Some(Internal(s,[a.[0];a.[1];Internal("",[a.[2];Leaf(leaf)])]))
                yield! 
                    enumerateAddToTree a.[0] leaf
                    |> Seq.choose id
                    |> Seq.map (fun t' -> Some(Internal(s,[t';a.[1];a.[2]])))
                yield!
                    enumerateAddToTree a.[1] leaf
                    |> Seq.choose id
                    |> Seq.map (fun t' -> Some(Internal(s,[a.[0];t';a.[2]])))
                yield!
                    enumerateAddToTree a.[2] leaf
                    |> Seq.choose id
                    |> Seq.map (fun t' -> Some(Internal(s,[a.[0];a.[1];t'])))
        | Leaf(l) -> yield None
        | Anon -> yield None
       }
       
let rec enumerateTrees (leaves:string[]) =
    if leaves.Length=3 then
        seq { yield Internal("",[Leaf(leaves.[0]);Leaf(leaves.[1]);Leaf(leaves.[2])]) }
    else 
        seq {
            for t in (enumerateTrees leaves.[1..]) do
                    yield! enumerateAddToTree t leaves.[0]
        }
        |> Seq.choose id

getData "eubt"
|> fun x -> x.Split(' ')
|> enumerateTrees
|> Seq.map printTree
|> Seq.iter (printfn "%s;")

// mend
type MendEst = { AA :float ; Aa: float ; aa: float }
let rec mend t =
    match t with 
    | Leaf(s) ->
        match s with 
        | "AA" -> { AA=1.; Aa=0.; aa=0. }
        | "Aa" -> { AA=0.; Aa=1.; aa=0. }
        | "aa" -> { AA=0.; Aa=0.; aa=1. }
    | Anon -> {AA=0.; Aa=0.; aa=0.}
    | Internal(s,l) ->
        l
        |> Array.ofList
        |> fun a -> 
            let m1=mend a.[0]
            let m2=mend a.[1]
            { 
              AA=m1.AA*m2.AA+m1.AA*0.5*m2.Aa+0.5*m1.Aa*m2.AA+0.25*m1.Aa*m2.Aa;
              Aa=m1.AA*0.5*m2.Aa+m1.AA*m2.aa+0.5*m1.Aa*m2.AA+0.5*m1.Aa*m2.Aa+0.5*m1.Aa*m2.aa+m1.aa*m2.AA+m1.aa*0.5*m2.Aa
              aa=0.25*m1.Aa*m2.Aa+0.5*m1.Aa*m2.aa+m1.aa*m2.aa+m1.aa*0.5*m2.Aa   
            }
            |> fun a -> 
                if abs (a.AA+a.Aa+a.aa-1.0)>1e-6 then eprintfn "Error: l:%A m1:%A m2:%A a:%A" l m1 m2 a;
                a


getData "mend"
|> splitNewline
|> Seq.take 1 |> Seq.exactlyOne
|> parseTree
|> mend
|> fun a -> printfn "%f %f %f" a.AA a.Aa a.aa


//chbp
// requires validation
// seems that we can't guarantee consistency.
// so maybe we need some ordering on splits in order to build the correct tree.
// we could enumerate all possible combinations and choose the first that validates- too expensive for >9 rows.
// http://evolution.berkeley.edu/evolibrary/article/0_0_0/phylogenetics_07
type Species = { Name:string; Characters:Map<int,int> }
let isSuperset s1 s2 = // is s2 a subset of s1?
    Seq.zip s1 s2
    |> Seq.fold (fun b (c1,c2)-> b&&(c1>=c2)) true
    
let complement s =
    s|>Seq.map (fun c -> 1-c)

let complementStr s =
    s|>Seq.map (fun c -> if c='1' then "0" else "1")|> String.concat ""

let toInts (s:string) = s|>Seq.map (fun c -> string c |> int)

let charToSplit taxa split =
    Seq.zip split taxa
    |> Seq.groupBy (fun (c,t)->c)
    |> Seq.map (fun (c,s)-> (int c,s|>Seq.map (fun (c,t)->t)|>Set.ofSeq))
    |> Map.ofSeq

[<StructuredFormatDisplay("{AsString}")>]
type phyloEdge = 
    { split:seq<int>;
      v1:int;
      v2:int }
    with 
    override this.ToString() = sprintf "(%d,%d,%s)" this.v1 this.v2 (this.split|>Seq.map string|>String.concat "")
    member m.AsString = m.ToString()

type phyloTree = int*list<phyloEdge> // vertex count and edges

let xorSplit s1 s2 =
        Seq.zip s1 s2
        |> Seq.map (fun (c1,c2) -> (c1+c2)%2)

// feels like we're on the right track.    
let rec addToTree (vcount,tree) _split =
    match tree with
    | [] -> 
        //eprintfn "%s" "Note: adding to empty tree"
        let e= { split=_split; v1=0; v2=1 }
        (2,[e])
    | head::tail ->
        //  is there a way to split head that's consistent with the new split?
        if isSuperset _split head.split then
            // add a new vertex vi splitting head, and another capturing the set not represented by the right end of head.
            let vi = vcount
            let vj = vi+1
            let e1 = { split=_split; v1= head.v1; v2=vi}
            let e2= { split = xorSplit _split head.split; v1=vi; v2=vj }
            let head' = {split=head.split; v1=vi; v2=head.v2 }
            (vcount+2,e1::e2::head'::tail)
        else if isSuperset head.split _split then
            let vi = vcount
            let vj = vi+1
            let e1 = { split = _split; v1=vi; v2=head.v2 }
            let e2= { split = xorSplit head.split _split; v1=vi; v2=vj }
            let head' = { split = head.split; v1=head.v1; v2=vi}
            (vcount+2,e1::e2::head'::tail)
        else
            let (v',t')=addToTree (vcount,tail) _split
            (v',head::t')

let toAdj tree =
    tree
    |> Seq.map (fun e -> [e.v1,e;e.v2,e])
    |> List.concat
    |> Seq.groupBy (fun (v,e)->v)
    |> Seq.map (fun (v,s)-> (v,Seq.map (fun (v,e)->e) s|>List.ofSeq))
    |> Map.ofSeq


// we have the edges of the tree and the adjacency list.  The tree is expressed in terms of splits.  Now we need to map the species to leaves of the tree.
let printAdj m =
    m
    |> Map.toSeq
    |> Seq.iter (fun (k,v)-> printfn "%d: %s" k (v|>Seq.map (fun e -> (if e.v1=k then e.v2 else e.v1)|>string)|>String.concat ","))

let validateTree (v,t) =
    toAdj t
    |> fun a -> printAdj a; a
    |> Map.toSeq
    |> Seq.map (fun (v,s)-> Seq.length s)
    |> Seq.countBy id
    |> eprintfn "Vertices: %d Counts: %A" v
    //t
    //|> fun (head::tail) -> (validateEdge tail head)
    t

// Buggy. Fix.
let assignTaxa taxa tree =
    let adj=toAdj tree
    let rec helper i parent e =
        let lr = e.split|>Seq.toArray|> fun a -> a.[i]
        let child = 
            if lr=0 then e.v1 else e.v2
        if child=parent then 
            eprintfn "%s" "Warning: child-parent conflict?"
            eprintfn "%A" e
            eprintfn "%A" adj.[child]
            eprintfn "%A" adj.[parent]
            parent
        else 
            let aa=adj.[child]
            //eprintfn "%A" aa
            if List.length aa = 1 then child
            else 
                // check each edge adjacent to child
                aa 
                |> List.map (helper i child)
                |> List.filter (fun v -> v<>child)
                |> fun (head::tail)->head
    taxa
    |> Seq.map (fun (t,i)-> (t,helper i -1 (List.head tree)))         

//#r @"C:\GitHub\automatic-graph-layout\GraphLayout\MSAGL\bin\Debug\Microsoft.Msagl.dll"
//#r @"C:\GitHub\automatic-graph-layout\GraphLayout\Drawing\bin\Debug\Microsoft.Msagl.Drawing.dll"
// #r @"System.Windows.Forms.dll"
// #r @"System.Windows.Forms.DataVisualization.dll"
// #r @"C:\GitHub\automatic-graph-layout\GraphLayout\tools\GraphViewerGDI\bin\Debug\Microsoft.Msagl.GraphViewerGdi.dll"
let render tree = 
    let form = new System.Windows.Forms.Form()
    //create a viewer object 
    let viewer = new Microsoft.Msagl.GraphViewerGdi.GViewer();
    //create a graph object 
    let graph = new Microsoft.Msagl.Drawing.PhyloTree() //new Microsoft.Msagl.Drawing.Graph("graph");
    //create the graph content 
    tree
    |> Seq.iter (fun e -> graph.AddEdge((string e.v1),(string e.v2))|>ignore)
   //bind the graph to the viewer 
    viewer.Graph <- graph;
    //associate the viewer with the form 
    form.SuspendLayout();
    viewer.Dock <- System.Windows.Forms.DockStyle.Fill;
    form.Controls.Add(viewer);
    form.ResumeLayout();
    //show the form 
    form.ShowDialog();



getData "chbp"
|> splitNewline
|> fun arr ->
    let taxa=arr.[0].Split(' ')|> Seq.mapi (fun i t -> (t,i))
    arr.[1..]
    |> Seq.map toInts
    |> Seq.fold addToTree (0,[])
   // |> fun t -> toAdj t, t
   // |> 
    |> validateTree
//    |> render
    
    |> assignTaxa taxa
    //|> Seq.iter (printfn "%A") *)

// TODO: fix this.
let rec makePhylogeny tree split  =
    match tree with
    | Anon -> Internal(split,[Leaf("");Leaf("")])
    | Internal(character,subtrees) ->
        let s1=toInts split
        let s2=toInts character
        let ss =
            (isSuperset s1 s2,
             isSuperset s1 (complement s2),
             isSuperset (complement s1) s2,
             isSuperset (complement s1) (complement s2),
             isSuperset s2 s1,
             isSuperset s2 (complement s1),
             isSuperset (complement s2) s1,
             isSuperset (complement s2) (complement s1))
        eprintfn "%s\n%s\n%A" split character ss // something wrong here.

        if isSuperset s1 s2 || isSuperset s1 (complement s2) then
            Internal(split,[Leaf("");Internal(character,subtrees)])
        else if isSuperset (complement s1) s2 || isSuperset (complement s1) (complement s2) then
            Internal(split,[Internal(character,subtrees);Leaf("")])
        else
            let a = subtrees|>Array.ofList
            //eprintfn "make: %d" a.Length
            if isSuperset s2 s1 || isSuperset s2 (complement s1) then // s2.[1] is superset of this character
                Internal(character,[a.[0]; makePhylogeny a.[1] split;] )
            else // must be true: isSuperset (complement s2) s1 || isSuperset (complement s2) (complement s1)
                eprintfn "isSuperset: %b" (isSuperset (complement s2) s1 || isSuperset (complement s2) (complement s1))
                Internal(character,[makePhylogeny a.[0] split; a.[1]] )
    | Leaf(l) -> 
        Internal(split,[Leaf("");Leaf("")])

let rec insertTaxa tree (t,i) =
    match tree with 
    | Anon -> Leaf(t)
    | Leaf(s) -> (Leaf(s+t+";"))   
    | Internal(character,subtrees) ->
        let c=character.[i]
        let a = subtrees|>Array.ofList
        //eprintfn "insert: %d" a.Length
        if c='0' then
            Internal(character,[insertTaxa a.[0] (t,i);a.[1]])
        else
            Internal(character,[a.[0];insertTaxa a.[1] (t,i)])    

let rec finalizeTree (t:NwckTree) : NwckTree option =
    match t with
    | Anon -> Some(Anon)
    | Leaf(s) -> 
        if (s="") then None
        else 
            s.Split([|';'|],StringSplitOptions.RemoveEmptyEntries)
            |> fun a ->             
                if (a.Length=1) then Some(Leaf(a.[0]))
                else Some(Internal("",[Leaf(a.[0]);seq{yield finalizeTree (Leaf(String.concat ";" a.[1..]))}|>Seq.choose id|>Seq.exactlyOne])) 
    | Internal(s,leaves) ->
        let result = leaves|>List.map finalizeTree|>List.choose id
        if (List.isEmpty result) then None
        else Some(Internal("",result))


getData "chbp" 
|> splitNewline
|> fun arr ->
    let taxa=arr.[0].Split(' ')|> Seq.mapi (fun i t -> (t,i))
    arr.[1..]
    |> List.ofArray
    |> List.fold makePhylogeny Anon
    |> fun t -> Seq.fold insertTaxa t taxa
    |> finalizeTree    
    |> fun t ->
        match t with 
        | None -> Anon
        | Some(s) -> s   
    |> printTree
    |> printfn "%s"

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

// sptd
let complementChar (s:string) =
    s
    |>Seq.map (fun c -> if c='0' then '1' else '0')
    |>Seq.map string
    |> String.concat ""

// a faster tree->character mapper
let rec getSubset =
        memoize (fun t ->
            match t with 
            | Leaf(s) -> s
            | Internal(s,l) ->   
                l
                |> List.map getSubset
                |> String.concat ",")

let taxaCount = memoize (fun mapping -> mapping|>Map.toSeq|>Seq.length)
    
let rec getCharacters (mapping:Map<string,int>) tree =
    
    let makeCharacter(s:string) =
        s.Split(',')
        |> Seq.map (fun t -> mapping.[t])
        |> Seq.sort
        |> Seq.map string
        |> String.concat ","

    match tree with
    | Leaf(s) ->  
        seq { yield makeCharacter s }
    | Internal(s,l) ->
        seq {
            for t in l do 
                yield makeCharacter (getSubset t)
                yield! getCharacters mapping t
        }

let getNonTrivialCharacters mapping tree =
    getCharacters mapping tree
    |> Seq.countBy id
    |> Seq.map (fun (k,c)->k)
    |> Seq.filter (fun s -> s.Split(',').Length>1)

let complementTaxa taxa (inSet:string) =
    let chkSet =
        inSet.Split(',')
        |> Set.ofArray
            
    taxa
    |> Map.toSeq
    |> Seq.map (fun (t,i)->i)
    |> Seq.filter (fun i-> not (Set.contains (string i) chkSet))
    |> Seq.sort
    |> Seq.map string
    |> String.concat ","

// functionally correct but too slow.
getData "sptd"
|> splitNewline
|> List.ofArray
|> fun (head::tail) ->
    let taxa = head.Split(' ')|> Seq.mapi (fun i t -> (t,i))|>Map.ofSeq
    let n = head.Split(' ').Length
    tail
    |> Seq.map parseTree
    |> Seq.map (getNonTrivialCharacters taxa)
    |> Array.ofSeq
    |> fun a ->
        let c0 = Set.ofSeq a.[0]
        a.[1]
        |> Seq.map (fun s -> if (Set.contains s c0) || (Set.contains (complementTaxa taxa s) c0) then 1 else 0)
        |> Seq.sum
        |> fun s-> 2*(n-3)-2*s
|> printfn "%A"

[<EntryPoint>]
let main argv = 
    printfn "%A" argv
    0 // return an integer exit code
