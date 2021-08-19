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

// eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
type State = any;

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  writableCategoriesEnabled:
    // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
    (state as any).config?.parameters?.annotations ?? false,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: (state as any).annoMatrix?.schema,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  userInfo: (state as any).userInfo,
}))
// eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
class Categories extends React.Component<{}, State> {
  // eslint-disable-next-line @typescript-eslint/ban-types --- FIXME: disabled temporarily on migrate to TS.
  constructor(props: {}) {
    super(props);
    this.state = {
      createAnnoModeActive: false,
      newCategoryText: "",
      categoryToDuplicate: null,
      expandedCats: new Set(),
    };
  }

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleCreateUserAnno = (e: any) => {
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleEnableAnnoMode = () => {
    this.setState({ createAnnoModeActive: true });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  handleDisableAnnoMode = () => {
    this.setState({
      createAnnoModeActive: false,
      categoryToDuplicate: null,
      newCategoryText: "",
    });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleModalDuplicateCategorySelection = (d: any) => {
    this.setState({ categoryToDuplicate: d });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  categoryNameError = (name: any) => {
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
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
    const { schema } = this.props;
    const allCategoryNames = schema.annotations.obs.columns.map(
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      (c: any) => c.name
    );
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleChange = (name: any) => {
    this.setState({ newCategoryText: name });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  handleSelect = (name: any) => {
    this.setState({ newCategoryText: name });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  instruction = (name: any) => {
    return labelPrompt(
      this.categoryNameError(name),
      "New, unique category name",
      ":"
    );
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types, @typescript-eslint/no-explicit-any -- - FIXME: disabled temporarily on migrate to TS.
  onExpansionChange = (catName: any) => {
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

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      createAnnoModeActive,
      categoryToDuplicate,
      newCategoryText,
      expandedCats,
    } = this.state;
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'writableCategoriesEnabled' does not exis... Remove this comment to see the full error message
    const { writableCategoriesEnabled, schema, userInfo } = this.props;
    /* all names, sorted in display order.  Will be rendered in this order */
    // @ts-expect-error ts-migrate(2554) FIXME: Expected 2 arguments, but got 1.
    const allCategoryNames = ControlsHelpers.selectableCategoryNames(
      schema
    ).sort();
    return (
      <div
        style={{
          padding: globals.leftSidebarSectionPadding,
        }}
      >
        <AnnoDialog
          // @ts-expect-error ts-migrate(2322) FIXME: Type '{ isActive: any; title: string; instruction:... Remove this comment to see the full error message
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
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ labelSuggestions: null; onChange: (name: a... Remove this comment to see the full error message
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
              // @ts-expect-error ts-migrate(2322) FIXME: Type '{ handleModalDuplicateCategorySelection: (d:... Remove this comment to see the full error message
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
        {/* eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS */}
        {allCategoryNames.map((catName: any) =>
          !schema.annotations.obsByName[catName].writable &&
          (schema.annotations.obsByName[catName].categories?.length > 1 ||
            !schema.annotations.obsByName[catName].categories) ? (
            <Category
              key={catName}
              // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
              metadataField={catName}
              onExpansionChange={this.onExpansionChange}
              isExpanded={expandedCats.has(catName)}
              createAnnoModeActive={createAnnoModeActive}
            />
          ) : null
        )}
        {/* WRITEABLE FIELDS */}
        {/* eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS */}
        {allCategoryNames.map((catName: any) =>
          schema.annotations.obsByName[catName].writable ? (
            <Category
              key={catName}
              // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
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
