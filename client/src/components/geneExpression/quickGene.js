import { H4, Icon, MenuItem } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { useState, useEffect, useRef, useMemo } from "react";
import fuzzysort from "fuzzysort";
import { Suggest } from "@blueprintjs/select";
import { useSelector, useDispatch } from "react-redux";

import Gene from "./gene";

import { postUserErrorToast } from "../framework/toasters";
import actions from "../../actions";

const usePrevious = (value) => {
  const ref = useRef();
  useEffect(() => {
    ref.current = value;
  });
  return ref.current;
};

function QuickGene() {
  const dispatch = useDispatch();

  const [isExpanded, setIsExpanded] = useState(true);
  const [geneNames, setGeneNames] = useState([]);
  const [, setStatus] = useState("pending");

  const { annoMatrix, userDefinedGenes, userDefinedGenesLoading } = useSelector(
    (state) => ({
        annoMatrix: state.annoMatrix,
        userDefinedGenes: state.controls.userDefinedGenes,
        userDefinedGenesLoading: state.controls.userDefinedGenesLoading,
      })
  );

  const prevProps = usePrevious({ annoMatrix });

  useEffect(async () => {
    if (!annoMatrix) return;
    if (annoMatrix !== prevProps?.annoMatrix) {
      const { schema } = annoMatrix;
      const varIndex = schema.annotations.var.index;

      setStatus("pending");
      try {
        const df = await annoMatrix.fetch("var", varIndex);
        setStatus("success");
        setGeneNames(df.col(varIndex).asArray());
      } catch (error) {
        setStatus("error");
        throw error;
      }
    }
  }, [annoMatrix, prevProps]);

  const handleExpand = () => setIsExpanded(!isExpanded);

  const renderGene = (fuzzySortResult, { handleClick, modifiers }) => {
    if (!modifiers.matchesPredicate) {
      return null;
    }
    /* the fuzzysort wraps the object with other properties, like a score */
    const geneName = fuzzySortResult.target;

    return (
      <MenuItem
        active={modifiers.active}
        disabled={modifiers.disabled}
        data-testid={`suggest-menu-item-${geneName}`}
        key={geneName}
        onClick={(g) =>
          /* this fires when user clicks a menu item */
          handleClick(g)
        }
        text={geneName}
      />
    );
  };

  const handleClick = (g) => {
    if (!g) return;
    const gene = g.target;
    if (userDefinedGenes.indexOf(gene) !== -1) {
      postUserErrorToast("That gene already exists");
    } else if (geneNames.indexOf(gene) === undefined) {
      postUserErrorToast("That doesn't appear to be a valid gene name.");
    } else {
      dispatch({ type: "single user defined gene start" });
      dispatch(actions.requestUserDefinedGene(gene));
      dispatch({ type: "single user defined gene complete" });
    }
  };

  const filterGenes = (query, genes) =>
    /* fires on load, once, and then for each character typed into the input */
    fuzzysort.go(query, genes, {
      limit: 5,
      threshold: -10000, // don't return bad results
    });

  const removeGene = (gene) => () => {
    dispatch({ type: "clear user defined gene", data: gene });
  };

  const QuickGenes = useMemo(() => userDefinedGenes.map((gene) => (
        <Gene
          key={`quick=${gene}`}
          gene={gene}
          removeGene={removeGene}
          quickGene
        />
      )), [userDefinedGenes]);

  return (
    <div style={{ width: "100%", marginBottom: "16px" }}>
      <H4
        role="menuitem"
        tabIndex="0"
        data-testclass="quickgene-heading-expand"
        onKeyPress={handleExpand}
        style={{
          cursor: "pointer",
        }}
        onClick={handleExpand}
      >
        Genes{" "}
        {isExpanded ? (
          <Icon icon={IconNames.CHEVRON_DOWN} />
        ) : (
          <Icon icon={IconNames.CHEVRON_RIGHT} />
        )}
      </H4>
      {isExpanded && (
        <>
          <div style={{ marginBottom: "8px" }}>
            <Suggest
              resetOnSelect
              closeOnSelect
              resetOnClose
              itemDisabled={userDefinedGenesLoading ? () => true : () => false}
              noResults={<MenuItem disabled text="No matching genes." />}
              onItemSelect={(g) => {
                /* this happens on 'enter' */
                handleClick(g);
              }}
              initialContent={<MenuItem disabled text="Enter a gene…" />}
              inputProps={{
                "data-testid": "gene-search",
                placeholder: "Quick Gene Search",
                leftIcon: IconNames.SEARCH,
                fill: true,
              }}
              inputValueRenderer={() => ""}
              itemListPredicate={filterGenes}
              itemRenderer={renderGene}
              items={geneNames || ["No genes"]}
              popoverProps={{ minimal: true }}
              fill
            />
          </div>
          {QuickGenes}
        </>
      )}
    </div>
  );
}

export default QuickGene;
