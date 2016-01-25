// Learn more about F# at http://fsharp.org
// See the 'F# Tutorial' project for more help.
open System
open System.IO

let (@@) folder filename = Path.Combine(folder,filename)
let data_root = @"c:\GitHub\rosalind\data"
let getData s = 
    File.ReadAllText(data_root@@(sprintf "rosalind_%s_1_dataset.txt" s)) |> fun x -> x.Trim()
let splitNewline (x:string) = x.Split([|'\r';'\n'|],StringSplitOptions.RemoveEmptyEntries)

type SuffixTree =    
    | Internal of Map<string,SuffixTree>
    | Leaf

let naiveMakeSuffixTree (s:string) =
    let matchStr (x:string) (y:string) =       
        // ugh slow
        let mutable i=0
        while i<x.Length && i<y.Length && x.[i]=y.[i] do i<-i+1
        x.Substring(0,i)
        (*
        Seq.zip x y
        |> Seq.takeWhile (fun (a,b)-> a=b)
        |> Seq.map (fun (a,_)-> string a)
        |> String.concat "" *)
        //|> fun sOut -> 
            //eprintfn "x:%s y:%s s:%s" x y sOut
        //    sOut
    let rec addToTree (t:SuffixTree) (sfx:string) : SuffixTree=
        match t with 
        | Leaf -> Map.empty|>Map.add sfx Leaf|>Internal
        | Internal(m) -> 
            m
            |> Map.toSeq
            |> Seq.map (fun (k,v) -> (k,matchStr k sfx))
            |> Seq.filter (fun (k,strmatch)-> strmatch.Length>0)
            |> fun s ->
                if Seq.isEmpty s then                
                    Map.add sfx Leaf m|>Internal
                else                
                    let (k,strMatch)=Seq.take 1 s|>Seq.exactlyOne
                    if strMatch.Length<k.Length then                        
                        let u=m
                        let v=m.[k]
                        let u'=Map.remove k m
                        let w=
                            Map.empty
                            |> Map.add (k.Substring(strMatch.Length)) v
                            |> Map.add (sfx.Substring(strMatch.Length)) Leaf
                            |> Internal
                        Map.add strMatch w u' |>Internal                       
                    else
                        Map.add k (addToTree m.[k] (sfx.Substring(strMatch.Length))) m|>Internal
    [0..s.Length-1]
    |> Seq.fold (fun tout i -> addToTree tout (s.Substring(i))) (Internal(Map.empty))


// max substrings of length k over |A|=a in a string of length n
let maxsubstrings (a:int) (n:int) (k:int) = 
    if k=1 then a
    else min (int (Math.Pow(float a,float k))) (n-k+1)
let m a n =
    [1..n]
    |> Seq.map (maxsubstrings a n)
    |> Seq.sum  

let sub s =
    let rec countSubs t =
        match t with
        | Leaf -> 0
        | Internal(m) ->
            m
            |> Map.toSeq
            |> Seq.map (fun (k,t')-> (k.Length)+(countSubs t')) 
            |> Seq.sum
    s // + "$" Is terminator necessary?  I don't think do since we're ok with prefixes getting mixed in. ??
    |> naiveMakeSuffixTree
    |> countSubs
      

let doLing s =
    (sub s, m 4 s.Length)
    |> fun (s,m)-> (float s)/(float m)
    

[<EntryPoint>]
let main argv = 
    
    getData "ling"
    |> splitNewline
    |> Seq.take 1 
    |> Seq.exactlyOne
    |> doLing
    |> printfn "%f"
    0 // return an integer exit code
