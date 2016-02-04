// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.
open System
open System.IO

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)


// this is a  much faster suffix tree that operates on s in-place, storing (ix,length) pairs as keys into the tree.  
// It's still not optimal but quite fast for these problem instances.
type SuffixTree2 =    
    | Internal2 of Map<int*int,SuffixTree2>
    | Leaf2

let naiveMakeSuffixTree2 (s:string) =
    let matchStr (a:int,alen:int) (b:int,blen:int) =
        let mutable i=0
        while i<alen && i<blen && s.[a+i]=s.[b+i] do
            i<-i+1
        i       
        (*
        Seq.zip x y
        |> Seq.takeWhile (fun (a,b)-> a=b)
        |> Seq.map (fun (a,_)-> string a)
        |> String.concat "" *)
        //|> fun sOut -> 
            //eprintfn "x:%s y:%s s:%s" x y sOut
        //    sOut
    let rec addToTree (t:SuffixTree2) (i:int,j:int) : SuffixTree2=
        match t with 
        | Leaf2 -> Map.empty|>Map.add (i,j) Leaf2|>Internal2
        | Internal2(m) -> 
            m
            |> Map.toSeq
            |> Seq.map (fun (k,v) -> (k,matchStr k (i,j)))
            |> Seq.filter (fun (k,strmatch)-> strmatch>0)
            |> fun s ->
                if Seq.isEmpty s then                
                    Map.add (i,j) Leaf2 m|>Internal2
                else                
                    let ((kix,klen),strMatch)=Seq.take 1 s|>Seq.exactlyOne
                    if strMatch<klen then                        
                        let u=m
                        let v=m.[kix,klen]
                        let u'=Map.remove (kix,klen) m
                        let w=
                            Map.empty
                            |> Map.add (kix+strMatch,klen-strMatch) v
                            |> Map.add (i+strMatch,j-strMatch) Leaf2
                            |> Internal2
                        Map.add (kix,strMatch) w u' |>Internal2                       
                    else
                        Map.add (kix,klen) (addToTree m.[kix,klen] (i+strMatch,j-strMatch)) m|>Internal2
    [0..s.Length-1]
    |> Seq.fold (fun tout i -> addToTree tout (i,s.Length-i)) (Internal2(Map.empty))
    |> fun t -> s,t

let rec enumerateSuffixEdges2 (s:string,t:SuffixTree2) =
    match t with
    | Leaf2 -> Seq.empty
    | Internal2 (m) ->
        seq {
            for ((i,j),v) in (Map.toSeq m) do
                yield! enumerateSuffixEdges2 (s,v)
                yield s.Substring(i,j)
        }
 
let printSuffixEdges2 t =
    enumerateSuffixEdges2 t
    |> Seq.iter (printfn "%s") 

// max substrings of length k over |A|=a in a string of length n
let maxsubstrings (a:int64) (n:int64) (k:int64) = 
    //eprintfn "%d %d %d" a n k
    if k=1L then a
    else if (float k)*Math.Log(float a)> Math.Log(float (n-k+1L)) then (n-k+1L) // what is a reasonable threshold?
    else
        // hacky way to find the smallest a^k' bigger than n-k+1L so we avoid overflow.  There's probably a faster/better way to do this.
        [1L..k]
        |> Seq.takeWhile (fun k -> Math.Pow(float a,float k)<float (n-k+1L)) 
        |> Seq.max
        |> fun k' ->
            min (int64 (Math.Pow(float a,float (k'+1L)))) (n-k+1L)
let m a n =    
    [1L..n]
    |> Seq.map (maxsubstrings a n)
    |> Seq.sum  

let sub s =
    let rec countSubs (s,t) =
        match t with
        | Leaf2 -> 0L
        | Internal2(m) ->
            m
            |> Map.toSeq
            |> Seq.map (fun ((kix,klen),t')-> (int64 klen)+(countSubs (s,t'))) 
            |> Seq.sum
    s // + "$" Is terminator necessary?  I don't think do since we're ok with prefixes getting mixed in. ??
    |> naiveMakeSuffixTree2
    |> countSubs
      

let doLing s =
    (sub s, m 4L (int64 s.Length))
    |> fun (s,m)-> (float s)/(float m)
    
// remove substrings
let filterStrings (ss:(int*string) seq) =
    ss
    |> Seq.sortBy (fun (_,s) -> s.Length)
    |> Array.ofSeq 
    |> fun a -> 
        [0..a.Length-1] 
        |> Seq.map (fun i ->
                [i+1..a.Length-1]
                |> Seq.fold (fun b j -> b||(snd a.[j]).Contains(snd a.[i])) false
                |> fun b -> if not b then Some(a.[i]) else None)
        |> Seq.choose id

let rec mrep (s:string,t:SuffixTree2) =
    let helper (kix,klen) m2 =
        mrep (s,m2)
        |> Seq.filter (fun (c,s')->c>1)
        |> Seq.map (fun (c,s')-> (c,s.Substring(kix,klen)+s'))
    
    match t with 
    | Leaf2 -> seq { yield (1,"") }
    | Internal2(m) ->
        m
        |> Map.toSeq
        |> Seq.map 
            (fun ((kix,klen),v) -> 
                match v with 
                | Leaf2 -> seq { yield (1,"") }
                | Internal2(m2)-> 
                    let r = 
                        helper (kix,klen) v
                      //  |> Seq.filter (fun (c,s)->c>1) unnecessary

                    if Seq.isEmpty r then // none of the children is a candidate to extend
                        seq {yield (Seq.length m2, s.Substring(kix,klen))}
                    else 
                        r // we can extend the string longer
                      //  |> Seq.map (fun (c,s')->c,s.Substring(kix,klen)+s') )
                      )
        |> Seq.concat
        |> filterStrings

let doMrep() =
    getData "mrep"
    |> splitNewline
    |> Seq.take 1 |> Seq.exactlyOne
    |> naiveMakeSuffixTree2
    |> mrep
    |> Seq.filter (fun (c,s) -> s.Length>=20)
    |> Seq.iter (fun (c,s)->printfn "%s" s)


[<EntryPoint>]
let main argv = 
    
    getData "ling"
    |> splitNewline
    |> Seq.take 1 
    |> Seq.exactlyOne
    |> doLing
    |> printfn "%f"
    0 // return an integer exit code
