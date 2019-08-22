import React from "react";
import _ from "lodash";
import { connect } from "react-redux";
import { FaChevronRight, FaChevronDown } from "react-icons/fa";
import { Button, Tooltip } from "@blueprintjs/core";

import * as globals from "../../globals";
import Value from "./value";
import sortedCategoryValues from "./util";

@connect(state => ({
  colorAccessor: state.colors.colorAccessor,
  categoricalSelection: state.categoricalSelection
}))
class Category extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      isChecked: true,
      isExpanded: false
    };
  }

  componentDidUpdate(prevProps) {
    const { categoricalSelection, metadataField } = this.props;
    if (categoricalSelection !== prevProps.categoricalSelection) {
      const cat = categoricalSelection[metadataField];
      const categoryCount = {
        // total number of categories in this dimension
        totalCatCount: cat.numCategoryValues,
        // number of selected options in this category
        selectedCatCount: _.reduce(
          cat.categoryValueSelected,
          (res, cond) => (cond ? res + 1 : res),
          0
        )
      };
      if (categoryCount.selectedCatCount === categoryCount.totalCatCount) {
        /* everything is on, so not indeterminate */
        this.checkbox.indeterminate = false;
        this.setState({ isChecked: true }); // eslint-disable-line react/no-did-update-set-state
      } else if (categoryCount.selectedCatCount === 0) {
        /* nothing is on, so no */
        this.checkbox.indeterminate = false;
        this.setState({ isChecked: false }); // eslint-disable-line react/no-did-update-set-state
      } else if (categoryCount.selectedCatCount < categoryCount.totalCatCount) {
        /* to be explicit... */
        this.checkbox.indeterminate = true;
        this.setState({ isChecked: false });
      }
    }
  }

  handleColorChange = () => {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "color by categorical metadata",
      colorAccessor: metadataField
    });
  };

  toggleAll() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter all of these",
      metadataField
    });
    this.setState({ isChecked: true });
  }

  toggleNone() {
    const { dispatch, metadataField } = this.props;
    dispatch({
      type: "categorical metadata filter none of these",
      metadataField
    });
    this.setState({ isChecked: false });
  }

  handleToggleAllClick() {
    const { isChecked } = this.state;
    // || this.checkbox.indeterminate === false
    if (isChecked) {
      this.toggleNone();
    } else if (!isChecked) {
      this.toggleAll();
    }
  }

  renderCategoryItems() {
    const { categoricalSelection, metadataField } = this.props;

    const cat = categoricalSelection[metadataField];
    const optTuples = sortedCategoryValues([...cat.categoryValueIndices]);
    return _.map(optTuples, (tuple, i) => (
      <Value
        optTuples={optTuples}
        key={tuple[1]}
        metadataField={metadataField}
        categoryIndex={tuple[1]}
        i={i}
      />
    ));
  }

  render() {
    const { isExpanded, isChecked } = this.state;
    const { metadataField, colorAccessor, categoricalSelection } = this.props;
    const { isTruncated } = categoricalSelection[metadataField];
    return (
      <div
        style={{
          maxWidth: globals.maxControlsWidth
        }}
        data-testclass="category"
        data-testid={`category-${metadataField}`}
      >
        <div
          style={{
            display: "flex",
            justifyContent: "space-between",
            alignItems: "baseline"
          }}
        >
          <div
            style={{
              display: "flex",
              justifyContent: "flex-start",
              alignItems: "baseline"
            }}
          >
            <label className="bp3-control bp3-checkbox">
              <input
                data-testclass="category-select"
                data-testid={`category-select-${metadataField}`}
                onChange={this.handleToggleAllClick.bind(this)}
                ref={el => {
                  this.checkbox = el;
                  return el;
                }}
                checked={isChecked}
                type="checkbox"
              />
              <span className="bp3-control-indicator" />
              {""}
            </label>

            <span
              data-testid={`category-expand-${metadataField}`}
              style={{
                cursor: "pointer",
                display: "inline-block"
              }}
              onClick={() => {
                this.setState({ isExpanded: !isExpanded });
              }}
            >
              {metadataField}
              {isExpanded ? (
                <FaChevronDown
                  data-testclass="category-expand-is-expanded"
                  style={{ fontSize: 10, marginLeft: 5 }}
                />
              ) : (
                <FaChevronRight
                  data-testclass="category-expand-is-not-expanded"
                  style={{ fontSize: 10, marginLeft: 5 }}
                />
              )}
            </span>
          </div>
          <Tooltip
            content="Use as color scale"
            position="bottom"
            hoverOpenDelay={globals.tooltipHoverOpenDelay}
          >
            <Button
              data-testclass="colorby"
              data-testid={`colorby-${metadataField}`}
              onClick={this.handleColorChange}
              active={colorAccessor === metadataField}
              intent={colorAccessor === metadataField ? "primary" : "none"}
              icon="tint"
            />
          </Tooltip>
        </div>
        <div style={{ marginLeft: 26 }}>
          {isExpanded ? this.renderCategoryItems() : null}
        </div>
        <div>
          {isExpanded && isTruncated ? (
            <p style={{ paddingLeft: 15 }}>... truncated list ...</p>
          ) : null}
        </div>
      </div>
    );
  }
}

export default Category;
