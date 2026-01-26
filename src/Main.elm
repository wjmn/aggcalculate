module Main exposing (..)

import Browser
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (..)
import List.Extra



---- MODEL ----


type alias PairedPeak =
    { aPeak : Float
    , tPeak : Float
    }


type EvaluatedResult
    = FailedNumberUniqueAlleles (List Float)
    | FailedAlleleSizeMatching (List Float)
    | PassedCheckWithScore Float


rank : EvaluatedResult -> Float
rank result =
    case result of
        FailedNumberUniqueAlleles _ ->
            10000.0

        FailedAlleleSizeMatching _ ->
            10000.0

        PassedCheckWithScore score ->
            score

type alias CalculatedResult =
    { pairedPeaks : List PairedPeak
    , pairedDistances : List ( Float, Float )
    , calculatedAlleleSizes : List Float
    , evaluation : EvaluatedResult
    }


type alias Model =
    { allele1Size : Int
    , allele2Size : Int
    , aBase : Int
    , tBase : Int
    , aPeaks : List Float
    , tPeaks : List Float
    , results : List CalculatedResult
    }


init : ( Model, Cmd Msg )
init =
    ( { allele1Size = 0
      , allele2Size = 0
      , aBase = 150
      , tBase = 151
      , aPeaks = []
      , tPeaks = []
      , results = []
      }
    , Cmd.none
    )



-- Produce a list of all possible permutations of a list


inserts : a -> List a -> List (List a)
inserts e xs =
    case xs of
        [] ->
            [ [ e ] ]

        head :: rest ->
            (e :: head :: rest) :: (inserts e rest |> List.map (\ys -> head :: ys))


permutations : List a -> List (List a)
permutations list =
    List.foldl (\x accum -> List.concatMap (inserts x) accum) [ [] ] list


-- If a list is smaller than size n, produce all possible lists of size n composed of the elements of the list
fillToSize : Int -> List a -> List (List a)
fillToSize n xs =
    let
        -- recursive helper function to fill the list to size n
        recHelper remainingIterations accumulatedLLists = 
            if remainingIterations <= 0 then
                accumulatedLLists
            else
                let
                    newLists =
                        List.concatMap (\l -> List.map (\x -> l ++ [ x ]) xs) accumulatedLLists
                in
                recHelper (remainingIterations - 1) newLists
    in
    
    if List.length xs >= n then
        [ xs ]
    else
        recHelper (n - List.length xs) [xs]

toGroupedString : Int -> List Int -> String
toGroupedString size interruptions = 
    case interruptions of 
        [] -> 
            "(CGG)" ++ String.fromInt size 
        head :: rest -> 
            let
                initials = 
                    head :: (List.map2 (-) (List.drop 1 interruptions) interruptions)
                    |> List.map (\x -> "(CGG)" ++ String.fromInt (x-1) ++ " AGG")
                    |> String.join " "
                remainder = 
                    " (CGG)" ++ String.fromInt (size - Maybe.withDefault 0 (List.maximum interruptions))
            in
            initials ++ remainder

---- UPDATE ----


type Msg
    = ChangedAllele1Size Int
    | ChangedAllele2Size Int
    | ChangedABase Int
    | ChangedTBase Int
    | ClickedAddAPeak
    | ChangedAPeak Int Float
    | RemovedAPeak Int
    | ClickedAddTPeak
    | ChangedTPeak Int Float
    | RemovedTPeak Int
    | NoOp


calculateFromModel : Model -> List CalculatedResult
calculateFromModel model =
    if model.aPeaks == [] && model.tPeaks == [] then
        []
    else
        let
            maxPeakCount =
                Basics.max (List.length model.aPeaks) (List.length model.tPeaks)
            filledAPeaks =
                fillToSize maxPeakCount model.aPeaks
            filledTPeaks =
                fillToSize maxPeakCount model.tPeaks
            -- Produce all possible pairs of A and T peaks
            aPeakPermutations = 
                List.concatMap permutations filledAPeaks
            pairedPeakList =
                filledTPeaks
                    |> List.concatMap (\tPeaks -> List.map (\aPerm -> List.map2 (\aPeak tPeak -> { aPeak = aPeak, tPeak = tPeak }) aPerm tPeaks) aPeakPermutations)

            evaluatePairedPeaks pairedPeaks =
                let
                    pairedDistances =
                        -- List.map (\p -> ( (p.aPeak - toFloat model.aBase) / 3, (p.tPeak - toFloat model.tBase) / 3 )) pairedPeaks
                        List.map (\p -> ( (p.aPeak - 164.2734995) / 2.872985178, (p.tPeak - 165.3278104) / 2.907717472 )) pairedPeaks -- empirical

                    calculatedAlleleSizes =
                        List.map (\( aDist, tDist ) -> aDist + tDist - 1) pairedDistances

                    sortedAlleles =
                        List.sort calculatedAlleleSizes

                    -- Find the number of gaps between consecutive alleles that exceed 3
                    gaps =
                        List.length <| List.filter (\( a, b ) -> b - a > 3) (List.map2 Tuple.pair sortedAlleles (List.drop 1 sortedAlleles))

                    -- For each calculated allele size, find the minimum of that size with each model allele size
                    minDiffs =
                        List.map (\size -> Basics.min (abs (toFloat model.allele1Size - size)) (abs (toFloat model.allele2Size - size))) calculatedAlleleSizes

                    evaluatedResult =
                        if gaps > 1 then
                            FailedNumberUniqueAlleles sortedAlleles

                        else if List.any (\diff -> diff > 2) minDiffs then
                            FailedAlleleSizeMatching sortedAlleles

                        else
                            PassedCheckWithScore (List.sum minDiffs / toFloat (List.length minDiffs))
                in
                { pairedPeaks = List.sortBy .aPeak pairedPeaks
                , pairedDistances = List.sort pairedDistances
                , calculatedAlleleSizes = List.sort calculatedAlleleSizes
                , evaluation = evaluatedResult
                }
        in
        pairedPeakList
            |> List.map evaluatePairedPeaks
            |> List.Extra.unique
            |> List.sortBy (.evaluation >> rank)

withNewCalculatedResult : Model -> Model
withNewCalculatedResult model =
    let
        newResult =
            calculateFromModel model
    in
    { model | results = newResult }


withCmd : Cmd Msg -> Model -> ( Model, Cmd Msg )
withCmd cmd model =
    ( model, cmd )


withCmdNone : Model -> ( Model, Cmd Msg )
withCmdNone model =
    ( model, Cmd.none )


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    case msg of
        ChangedAllele1Size size ->
            { model | allele1Size = size }
                |> withNewCalculatedResult
                |> withCmdNone

        ChangedAllele2Size size ->
            { model | allele2Size = size }
                |> withNewCalculatedResult
                |> withCmdNone

        ChangedABase size ->
            { model | aBase = size }
                |> withNewCalculatedResult
                |> withCmdNone

        ChangedTBase size ->
            { model | tBase = size }
                |> withNewCalculatedResult
                |> withCmdNone

        ClickedAddAPeak ->
            { model | aPeaks = model.aPeaks ++ [ 0 ] }
                |> withCmdNone

        ChangedAPeak index value ->
            let
                newPeaks =
                    List.indexedMap
                        (\i peak ->
                            if i == index then
                                value

                            else
                                peak
                        )
                        model.aPeaks
            in
            { model | aPeaks = newPeaks }
                |> withNewCalculatedResult
                |> withCmdNone

        RemovedAPeak index ->
            let
                newPeaks =
                    List.indexedMap Tuple.pair model.aPeaks
                        |> List.filterMap
                            (\( i, peak ) ->
                                if i /= index then
                                    Just peak

                                else
                                    Nothing
                            )
            in
            { model | aPeaks = newPeaks }
                |> withNewCalculatedResult
                |> withCmdNone

        ClickedAddTPeak ->
            { model | tPeaks = model.tPeaks ++ [ 0 ] }
                |> withCmdNone

        ChangedTPeak index value ->
            let
                newPeaks =
                    List.indexedMap
                        (\i peak ->
                            if i == index then
                                value

                            else
                                peak
                        )
                        model.tPeaks
            in
            { model | tPeaks = newPeaks }
                |> withNewCalculatedResult
                |> withCmdNone

        RemovedTPeak index ->
            let
                newPeaks =
                    List.indexedMap Tuple.pair model.tPeaks
                        |> List.filterMap
                            (\( i, peak ) ->
                                if i /= index then
                                    Just peak

                                else
                                    Nothing
                            )
            in
            { model | tPeaks = newPeaks }
                |> withNewCalculatedResult
                |> withCmdNone

        NoOp ->
            model
                |> withCmdNone



---- VIEW ----


view : Model -> Html Msg
view model =

    let
        (bestResult, restOfResults) = 
            case model.results of 
                r :: tail-> 
                    case r.evaluation of 
                        PassedCheckWithScore result -> 
                            (div [class "main-success"] [viewResult model True r], tail)
                        _ -> 
                            (div [class "main-failed"] [text "The input parameters do not correspond with a valid AGG result. See below for all possibilities excluded during this calculation for troubleshooting."], model.results)
                _ -> 
                    (div [class "main-failed"] [text "The input parameters do not correspond with a valid AGG result. See below for all possibilities excluded during this calculation for troubleshooting."], model.results)

    in
    div [ class "container" ]
        [ div [ class "left-panel" ]
            [ h1 [class "title"] [ text "AGGCalculate v1.0" ]
            , h2 [class "subtitle"] [text "Input the repeat sizes of both alleles:"]
            , div [class "input-group"]
                [ label [] [ text "Allele 1 Size (rpts)" ]
                , input [ Html.Attributes.type_ "number", Html.Attributes.value (String.fromInt model.allele1Size), onInput (ChangedAllele1Size << Maybe.withDefault 0 << String.toInt) ] []
                ]
            , div [class "input-group"]
                [ label [] [ text "Allele 2 Size (rpts)" ]
                , input [ Html.Attributes.type_ "number", Html.Attributes.value (String.fromInt model.allele2Size), onInput (ChangedAllele2Size << Maybe.withDefault 0 << String.toInt) ] []
                ]
            -- , div [ class "input-group" ]
            --     [ label [] [ text "A Base Size (bp)" ]
            --     , input [ Html.Attributes.type_ "number", Html.Attributes.value (String.fromInt model.aBase), onInput (ChangedABase << Maybe.withDefault 0 << String.toInt) ] []
            --     ]
            -- , div [ class "input-group" ]
            --     [ label [] [ text "T Base Size (bp)" ]
            --     , input [ Html.Attributes.type_ "number", Html.Attributes.value (String.fromInt model.tBase), onInput (ChangedTBase << Maybe.withDefault 0 << String.toInt) ] []
            --     ]
            , br [] []
            , p [] [text "Input the peaks detected in the A-primed assay:"]
            , div [class "input-group"]
                [ label [] [text "A-Primed PCR Peaks (bp)"]
                , div [style "width" "100%"]
                    (List.indexedMap
                        (\index peak ->
                            div [class "row-item"]
                                [ input [ Html.Attributes.type_ "number", Html.Attributes.value (String.fromFloat peak), onInput (ChangedAPeak index << Maybe.withDefault 0 << String.toFloat) ] []
                                , button [ onClick (RemovedAPeak index) ] [ text "Remove" ]
                                ]
                        )
                        model.aPeaks
                    )
                , Html.button [ onClick ClickedAddAPeak ] [ text "Add A-Primer Reaction Peak" ]
                ]
            , br [] []
            , p [] [text "Input the peaks detected in the T-primed assay:"]
            , div [class "input-group"]
                [ label [] [text "T-Primed PCR Peaks (bp)"]
                , div [style "width" "100%"]
                    (List.indexedMap
                        (\index peak ->
                            div [class "row-item"]
                                [ input [ Html.Attributes.type_ "number", Html.Attributes.value (String.fromFloat peak), onInput (ChangedTPeak index << Maybe.withDefault 0 << String.toFloat) ] []
                                , button [ onClick (RemovedTPeak index) ] [ text "Remove" ]
                                ]
                        )
                        model.tPeaks
                    )
                , Html.button [ onClick ClickedAddTPeak ] [ text "Add T-Primer Reaction Peak" ]
                ]
            ]
        , div [ class "right-panel" ]
            [ div []
                [ h1 [class "results-title"] [text "Results" ]
                , div [] [bestResult]
                , Html.ul [] (List.map (viewResult model False) restOfResults)
                ]
            ]
        ]


viewResult : Model -> Bool -> CalculatedResult -> Html Msg
viewResult model isMainResult result =
    case result.evaluation of
        FailedNumberUniqueAlleles alleles ->
            details [class "failed"]
                [ summary [] [ text "Possibility Excluded: Number of unique alleles is too high." ]
                , div [] 
                    [ text "Paired peaks: ", 
                      text <| String.join ", " (List.map (\p -> "(" ++ String.fromFloat p.aPeak ++ ", " ++ String.fromFloat p.tPeak ++ ")") result.pairedPeaks)
                    , br [] []
                    , text "Paired Distances: ",
                      text <| String.join ", " (List.map (\(a, t) -> "(" ++ String.fromFloat a ++ ", " ++ String.fromFloat t ++ ")") result.pairedDistances)
                    , br [] []
                    , text "Calculated Allele Sizes: ",
                      text <| String.join ", " (List.map String.fromFloat result.calculatedAlleleSizes)
                    , br [] []
                    ]
                ]

        FailedAlleleSizeMatching alleles ->
            details [class "failed"]
                [ summary [] [ text "Possibility Excluded: Allele sizes do not match." ]
                , div [] 
                    [ text "Paired peaks: ", 
                      text <| String.join ", " (List.map (\p -> "(" ++ String.fromFloat p.aPeak ++ ", " ++ String.fromFloat p.tPeak ++ ")") result.pairedPeaks)
                    , br [] []
                    , text "Paired Distances: ",
                      text <| String.join ", " (List.map (\(a, t) -> "(" ++ String.fromFloat a ++ ", " ++ String.fromFloat t ++ ")") result.pairedDistances)
                    , br [] []
                    , text "Calculated Allele Sizes: ",
                      text <| String.join ", " (List.map String.fromFloat result.calculatedAlleleSizes)
                    , br [] []
                    ]
                ]

        PassedCheckWithScore score ->

            let
                groupedPairedDistances = 
                    result.pairedDistances
                        |> List.map (\ (a, t) -> if abs (a + t - toFloat model.allele1Size) > abs (a + t - toFloat model.allele2Size) then (model.allele2Size, round t) else (model.allele1Size, round t))
                        |> List.append [(model.allele1Size, 0), (model.allele2Size, 0)]
                        |> List.Extra.gatherEqualsBy Tuple.first
                        |> List.map (\( el1, pairs ) -> ( Tuple.first el1, List.map Tuple.second (el1 :: pairs) |> List.filter (\x -> x > 0) |> List.sort))
                viewSizeAndInterrupts ( size, interrupts) = 
                    let
                        numInterrupts = List.length interrupts
                    in
                    div []
                        [ text <| "The ", span [class "bold"] [text <| String.fromInt size ++ " repeat"], text <| " allele has ", span [class "bold"] [text <| String.fromInt numInterrupts ++ " AGG interruption(s)"], text <| " at position(s): " ++ 
                            String.join ", " (List.map String.fromInt interrupts)
                        ]
                viewSizeAndInterruptsSummary (size, interrupts) = 
                    let
                        numInterrupts = List.length interrupts
                    in
                    div []
                        [ text <| "The ", span [class "bold"] [text <| String.fromInt size ++ " repeat"], text <| " allele has ", span [class "bold"] [text <| String.fromInt numInterrupts ++ " AGG interruption(s)"], text " = ",
                            span [class "bold"] [text (toGroupedString size interrupts)]
                        ]
                headerText = 
                    if isMainResult then 
                        "Most Likely Valid AGG Result"
                    else
                        "Other Possible Valid AGG Result"

            in
            
            div [class "passed", classList [("main-result", isMainResult)] ]
                [ div [class "passed-header", classList [("main-result", isMainResult)]] [ div [class "passed-header-text"] [text headerText ]]
                , div [class "passed-body", classList[("main-result", isMainResult)]]
                [div [] (List.map viewSizeAndInterruptsSummary groupedPairedDistances)
                    , details  [class "passed-details"] 
                        [summary [class "passed-details-summary"] [text "Details"]
                        , div [class "passed-details-additional"] 
                        [ text <| "Evaluation Score (lower is better): " ++ String.fromFloat score
                        , br [] []
                        , text "Paired peaks: ", 
                        text <| String.join ", " (List.map (\p -> "(" ++ String.fromFloat p.aPeak ++ ", " ++ String.fromFloat p.tPeak ++ ")") result.pairedPeaks)
                        , br [] []
                        , text "Paired Distances: ",
                        text <| String.join ", " (List.map (\(a, t) -> "(" ++ String.fromFloat a ++ ", " ++ String.fromFloat t ++ ")") result.pairedDistances)
                        , br [] []
                        , text "Calculated Allele Sizes: ",
                        text <| String.join ", " (List.map String.fromFloat result.calculatedAlleleSizes)
                        , br [] []
                        , div [] (List.map viewSizeAndInterrupts groupedPairedDistances)
                        ]
                    ]
                ]
                ]



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.element
        { view = view
        , init = \_ -> init
        , update = update
        , subscriptions = always Sub.none
        }
