import React from "react";
import { AnchorButton, Tooltip, Position } from "@blueprintjs/core";
import { connect } from "react-redux";
import * as globals from "../../globals";
import Category from "./category";
import { AnnotationsHelpers, ControlsHelpers } from "../../util/stateManager";
import AnnoDialog from "../annoDialog";
import AnnoSelect from "./annoSelect";
import LabelInput from "../labelInput";
import { labelPrompt } from "./labelUtil";
import actions from "../../actions";

@connect((state) => ({
  writableCategoriesEnabled: state.config?.parameters?.annotations ?? false,
  schema: state.annoMatrix?.schema,
  userInfo: state.userInfo,
}))
class Categories extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      createAnnoModeActive: false,
      newCategoryText: "",
      categoryToDuplicate: null,
      expandedCats: new Set(),
    };
  }

  handleCreateUserAnno = (e) => {
    const { dispatch } = this.props;
    const { newCategoryText, categoryToDuplicate } = this.state;
    dispatch(
      actions.annotationCreateCategoryAction(
        newCategoryText,
        categoryToDuplicate
      )
    );
    this.setState({
      createAnnoModeActive: false,
      categoryToDuplicate: null,
      newCategoryText: "",
    });
    e.preventDefault();
  };

  handleEnableAnnoMode = () => {
    this.setState({ createAnnoModeActive: true });
  };

  handleDisableAnnoMode = () => {
    this.setState({
      createAnnoModeActive: false,
      categoryToDuplicate: null,
      newCategoryText: "",
    });
  };

  handleModalDuplicateCategorySelection = (d) => {
    this.setState({ categoryToDuplicate: d });
  };

  categoryNameError = (name) => {
    /*
    return false if this is a LEGAL/acceptable category name or NULL/empty string,
    or return an error type.
    */

    /* allow empty string */
    if (name === "") return false;

    /*
    test for uniqueness against *all* annotation names, not just the subset
    we render as categorical.
    */
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map((c) => c.name);

    /* check category name syntax */
    const error = AnnotationsHelpers.annotationNameIsErroneous(name);
    if (error) {
      return error;
    }

    /* disallow duplicates */
    if (allCategoryNames.indexOf(name) !== -1) {
      return "duplicate";
    }

    /* otherwise, no error */
    return false;
  };

  handleChange = (name) => {
    this.setState({ newCategoryText: name });
  };

  handleSelect = (name) => {
    this.setState({ newCategoryText: name });
  };

  instruction = (name) => labelPrompt(
      this.categoryNameError(name),
      "New, unique category name",
      ":"
    );

  onExpansionChange = (catName) => {
    const { expandedCats } = this.state;
    if (expandedCats.has(catName)) {
      const _expandedCats = new Set(expandedCats);
      _expandedCats.delete(catName);
      this.setState({ expandedCats: _expandedCats });
    } else {
      const _expandedCats = new Set(expandedCats);
      _expandedCats.add(catName);
      this.setState({ expandedCats: _expandedCats });
    }
  };

  render() {
    const {
      createAnnoModeActive,
      categoryToDuplicate,
      newCategoryText,
      expandedCats,
    } = this.state;
    const { writableCategoriesEnabled, schema, userInfo } = this.props;
    /* all names, sorted in display order.  Will be rendered in this order */
    const allCategoryNames =
      ControlsHelpers.selectableCategoryNames(schema).sort();

    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding,
        }}
      >
        <AnnoDialog
          isActive={createAnnoModeActive}
          title="Create new category"
          instruction={this.instruction(newCategoryText)}
          cancelTooltipContent="Close this dialog without creating a category."
          primaryButtonText="Create new category"
          primaryButtonProps={{ "data-testid": "submit-category" }}
          text={newCategoryText}
          validationError={this.categoryNameError(newCategoryText)}
          handleSubmit={this.handleCreateUserAnno}
          handleCancel={this.handleDisableAnnoMode}
          annoInput={
            <LabelInput
              labelSuggestions={null}
              onChange={this.handleChange}
              onSelect={this.handleSelect}
              inputProps={{
                "data-testid": "new-category-name",
                leftIcon: "tag",
                intent: "none",
                autoFocus: true,
              }}
              newLabelMessage="New category"
            />
          }
          annoSelect={
            <AnnoSelect
              handleModalDuplicateCategorySelection={
                this.handleModalDuplicateCategorySelection
              }
              categoryToDuplicate={categoryToDuplicate}
              allCategoryNames={allCategoryNames}
            />
          }
        />

        {writableCategoriesEnabled ? (
          <div style={{ marginBottom: 10 }}>
            <Tooltip
              content={
                userInfo.is_authenticated
                  ? "Create a new category"
                  : "You must be logged in to create new categorical fields"
              }
              position={Position.RIGHT}
              boundary="viewport"
              hoverOpenDelay={globals.tooltipHoverOpenDelay}
              modifiers={{
                preventOverflow: { enabled: false },
                hide: { enabled: false },
              }}
            >
              <AnchorButton
                type="button"
                data-testid="open-annotation-dialog"
                onClick={this.handleEnableAnnoMode}
                intent="primary"
                disabled={!userInfo.is_authenticated}
              >
                Create new <strong>category</strong>
              </AnchorButton>
            </Tooltip>
          </div>
        ) : null}

        {/* READ ONLY CATEGORICAL FIELDS */}
        {/* this is duplicative but flat, could be abstracted */}
        {allCategoryNames.map((catName) =>
          !schema.annotations.obsByName[catName].writable &&
          (schema.annotations.obsByName[catName].categories?.length > 1 ||
            !schema.annotations.obsByName[catName].categories) ? (
            <Category
              key={catName}
              metadataField={catName}
              onExpansionChange={this.onExpansionChange}
              isExpanded={expandedCats.has(catName)}
              createAnnoModeActive={createAnnoModeActive}
            />
          ) : null
        )}
        {/* WRITEABLE FIELDS */}
        {allCategoryNames.map((catName) =>
          schema.annotations.obsByName[catName].writable ? (
            <Category
              key={catName}
              metadataField={catName}
              onExpansionChange={this.onExpansionChange}
              isExpanded={expandedCats.has(catName)}
              createAnnoModeActive={createAnnoModeActive}
            />
          ) : null
        )}
      </div>
    );
  }
}

export default Categories;
